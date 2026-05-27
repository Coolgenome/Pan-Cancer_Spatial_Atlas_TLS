import sys
import PIL
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import scipy.ndimage as ndi
import pickle
from collections import Counter
import cv2
import csv
import os
import pandas as pd

prefix = sys.argv[1]
mask_suffix = sys.argv[2]
out_suffix = sys.argv[3]
csv_path = sys.argv[4]
mask_path= prefix + mask_suffix

out_path = prefix + out_suffix
os.makedirs(out_path, exist_ok=True)

PIL.Image.MAX_IMAGE_PIXELS = None

# Load your image and mask (assuming prefix and mask_path are defined)
he_image = Image.open(prefix+'he.jpg')
mask_tls = Image.open(mask_path).convert('L')  # Convert mask to grayscale

# Resize mask to match image dimensions
mask_tls_resized = mask_tls.resize(he_image.size)

he_image_np = np.array(he_image)
mask_tls_resized_np = np.array(mask_tls_resized)

mask_tls_resized_np = mask_tls_resized_np > 0
he_image_np[~mask_tls_resized_np] = 0

mask_tls_array = np.array(mask_tls_resized)
labeled_array, num_features = ndi.label(mask_tls_array > 0)

df = pd.read_csv(csv_path)
tls_ids = set(df['TLS_ID'])
# Load cluster data from the pickle file
#may modify cluster_pickle path if necessary
cluster_pickle_path = prefix + '/tls_cluster/nclusters005/labels.pickle'
with open(cluster_pickle_path, 'rb') as f:
    cluster_data = pickle.load(f)

# Filter out TLS_IDs that correspond to 'TO-P'
top_ids = set(df[df['Manual'] == 'TO-P']['TLS_ID'])
top_labeled_indices = [label for label in top_ids if label in np.unique(labeled_array)]
top_labeled_indices = sorted(top_labeled_indices)

resized_cluster_data = cv2.resize(cluster_data, (labeled_array.shape[1], labeled_array.shape[0]), interpolation=cv2.INTER_NEAREST)

# Initialize a dictionary to store positions for each TO-P TLS
top_tls_positions = {tls_id: [] for tls_id in top_labeled_indices}

# Find positions in labeled_array corresponding to each TO-P TLS
for tls_id in top_labeled_indices:
    positions = np.argwhere(labeled_array == tls_id)
    top_tls_positions[tls_id].extend(positions)

# Extract corresponding cluster information
top_tls_resized_cluster_data = {}
for tls_id, positions in top_tls_positions.items():
    clusters = [resized_cluster_data[pos[0], pos[1]] for pos in positions]
    fifth_cluster_positions = [pos for pos, cluster in zip(positions, clusters) if cluster == 4]
    top_tls_resized_cluster_data[tls_id] = fifth_cluster_positions


# Create a copy of labeled_array to modify
modified_labeled_array = labeled_array.copy()

# Iterate through each TO-P TLS ID and its positions
for tls_id, positions in top_tls_positions.items():
    for pos in positions:
        x, y = pos
        # Use resized_labeled_array for checking positions
        if labeled_array[x, y] == tls_id and resized_cluster_data[x, y] != 4:
            # Set to 0 in modified_labeled_array
            modified_labeled_array[x, y] = 0

# Create a mask from the modified labeled array
# The mask is True where modified_labeled_array is not zero
modified_mask = modified_labeled_array > 0
modified_he_image_np = np.where(modified_mask[..., np.newaxis], he_image_np, 0)

def mkdir(path):
    dirname = os.path.dirname(path)
    if dirname != '':
        os.makedirs(dirname, exist_ok=True)

def save_image(img, filename):
    mkdir(filename)
    Image.fromarray(img).save(filename)
    print(filename)


# Initialize a dictionary to store the new labeled arrays and their features
new_labeled_arrays = {}

# Iterate over top_labeled_indices
for tls_id in top_labeled_indices:
    tls_mask = modified_labeled_array == tls_id
    tls_labeled_array, tls_num_features = ndi.label(tls_mask)
    new_labeled_arrays[tls_id] = {
        'labeled_array': tls_labeled_array,
        'num_features': tls_num_features}


# Create a new array based on the modified_labeled_array
combined_labeled_array = modified_labeled_array.copy()

# Initialize a variable to keep track of the highest label number used
max_label = combined_labeled_array.max()

# Iterate through each TLS ID and its labeled array in new_labeled_arrays
for tls_id, data in new_labeled_arrays.items():
    tls_labeled_array = data['labeled_array']
    tls_num_features = data['num_features']
    for i in range(1, tls_num_features + 1):
        max_label += 1
        segment_positions = np.where(tls_labeled_array == i)
        combined_labeled_array[segment_positions] = max_label

# Update the total number of features in the combined labeled array
combined_num_features = max_label

def save_pickle(x, filename):
    mkdir(filename)
    with open(filename, 'wb') as file:
        pickle.dump(x, file)
    print(filename)

combined_labeled_array_pickle_path = out_path+'combined_labeled_array.pickle'
save_pickle(combined_labeled_array, combined_labeled_array_pickle_path)
