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

sig_suffix = "markers/phenotype/raw/"
gene_suffix = "cnts-super/"

# Concatenate to create the full paths
mask_path= prefix + mask_suffix
sig_pickle_path = prefix + sig_suffix
gene_pickle_path = prefix + gene_suffix
out_path = prefix + out_suffix
os.makedirs(out_path, exist_ok=True)

plot_path = prefix + 'cnts-super-plots/'

PIL.Image.MAX_IMAGE_PIXELS = None

df = pd.read_csv(csv_path)

# Load your image and mask (assuming prefix and mask_path are defined)
he_image = Image.open(prefix+'he.jpg')
he_image_np = np.array(he_image)
mask_tls = Image.open(mask_path).convert('L')  # Convert mask to grayscale
mask_tls_resized = mask_tls.resize(he_image.size)
mask_tls_array = np.array(mask_tls_resized) > 0

# Label the mask
labeled_array, num_features = ndi.label(mask_tls_array)

# Load cluster data from the pickle file
#may modify cluster_pickle path if necessary
cluster_pickle_path = prefix + '/tls_cluster/nclusters005/labels.pickle'
with open(cluster_pickle_path, 'rb') as f:
    cluster_data = pickle.load(f)

manual_col = df['Manual'].astype(str).str.strip().str.upper()
tls_ids_df = set(df.loc[manual_col == 'TLS', 'TLS_ID'].astype(int))
tls_ids_in_original = sorted(tls_ids_df.intersection(np.unique(labeled_array)))

# Filter out TLS_IDs that correspond to 'TO-P'
top_ids = set(df[df['Manual'] == 'TO-P']['TLS_ID'])
# Find the labeled indices for 'TO-P'
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


def mkdir(path):
    dirname = os.path.dirname(path)
    if dirname != '':
        os.makedirs(dirname, exist_ok=True)
        
def save_image(img, filename):
    mkdir(filename)
    Image.fromarray(img).save(filename)
    print(filename)

# Create a mask from the modified labeled array
# The mask is True where modified_labeled_array is not zero
modified_mask = modified_labeled_array > 0
modified_he_image_np = np.where(modified_mask[..., np.newaxis], he_image_np, 0)

# Initialize a dictionary to store the new labeled arrays and their features
new_labeled_arrays = {}

# Iterate over top_labeled_indices
for tls_id in top_labeled_indices:
    # Create a new mask for the current TLS ID using the modified_labeled_array
    tls_mask = modified_labeled_array == tls_id

    # Label this new mask
    tls_labeled_array, tls_num_features = ndi.label(tls_mask)

    # Store the labeled array and the number of features
    new_labeled_arrays[tls_id] = {
        'labeled_array': tls_labeled_array,
        'num_features': tls_num_features
    }

# Create a new array based on the modified_labeled_array
combined_labeled_array = modified_labeled_array.copy()

# Initialize a variable to keep track of the highest label number used
max_label = combined_labeled_array.max()

# Iterate through each TLS ID and its labeled array in new_labeled_arrays
for tls_id, data in new_labeled_arrays.items():
    tls_labeled_array = data['labeled_array']
    tls_num_features = data['num_features']

    # Loop through each feature/segment for the current TLS ID
    for i in range(1, tls_num_features + 1):
        # Increment the max_label for a new unique label
        max_label += 1

        # Find positions where the tls_labeled_array has the label i
        segment_positions = np.where(tls_labeled_array == i)

        # Update these positions in combined_labeled_array with the new label
        combined_labeled_array[segment_positions] = max_label

# Update the total number of features in the combined labeled array
combined_num_features = max_label
tls_ids_in_combined = sorted(tls_ids_df.intersection(np.unique(combined_labeled_array)))

def segment_image(img,output_filename,label_ids,combined_labeled_array):
    segmented_images = []
    # Loop through each feature/segment for the current TLS ID
    for i in label_ids:
        # Create a mask for the current segment
        segment_mask = (combined_labeled_array == i).astype(np.uint8) * 255

        # Convert the NumPy array back to PIL image for cropping
        segment_mask_pil = Image.fromarray(segment_mask)

        # Apply this mask to the image
        img.putalpha(segment_mask_pil)

        # Get the bounding box and crop the image
        bbox = segment_mask_pil.getbbox()
        if bbox:
            cropped_segment = img.crop(bbox)

            # Store or process the cropped segment
            segmented_images.append(cropped_segment)

            # Save each segment with a unique filename including the TLS ID and segment index
            filename = f'{out_path}{output_filename}_segmented_image_{i}.png'
            cropped_segment.save(f'{filename}')
    return segmented_images

genes = ['CD3D','MS4A1','CR2', 'FCER2', 'BCL6', 'MKI67']

plot_files = [f for f in os.listdir(plot_path) if f.endswith('.png')]

existing_png_files = []

# Check if each gene has a corresponding pickle file and load the data
for gene in genes:
    plot_filename_gene = gene + '.png'
    if plot_filename_gene in plot_files:
        full_path_plot = os.path.join(plot_path, plot_filename_gene)
        # Open the image using PIL
        gene_he = Image.open(full_path_plot)
        # Resize the image to match 'he_image.size'
        gene_he_resized = gene_he.resize(he_image.size)
        segmented_images_TOP = segment_image(gene_he_resized, gene, tls_ids_in_combined, combined_labeled_array)
        segmented_images = segment_image(gene_he_resized, gene,tls_ids_in_original,labeled_array)

segmented_he_TOP = segment_image(he_image, "HE", tls_ids_in_combined, combined_labeled_array)
segmented_he = segment_image(he_image, "HE",tls_ids_in_original,labeled_array)

