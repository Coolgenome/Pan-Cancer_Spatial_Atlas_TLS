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


###number of pixels
segment_pixel_counts = {}

for i in range(1, combined_num_features + 1):
    # Create a mask for the current segment
    segment_mask = (combined_labeled_array == i).astype(np.uint8) * 255

    # Calculate the number of pixels in the current segment
    num_pixels = np.count_nonzero(segment_mask)

    # Store the pixel count in the dictionary
    segment_pixel_counts[i] = num_pixels

csv_pixel_filename = f'{out_path}segment_pixel_counts.csv'

with open(csv_pixel_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Writing the header
    writer.writerow(['Label', 'Pixel Count', 'Manual'])

    # Writing the data
    for segment_label, pixel_count in segment_pixel_counts.items():
         # Write only if pixel count is non-zero
        if pixel_count > 0:
            manual_value = 4 if segment_label > num_features else ''
            writer.writerow([segment_label, pixel_count, manual_value])

pickle_files = [
    f'{sig_pickle_path}tls.pickle',
    f'{sig_pickle_path}lymphoid.pickle'
]
# Dictionary to store the loaded data from each file
loaded_data = {}

# Iterate over the list of file paths and load each file
for file_path in pickle_files:
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
        loaded_data[file_path] = data

# Dictionary to store score records for each file
all_score_records = {}
all_score_records2 = {}
# Process each set of scores
for file_path, score in loaded_data.items():
    score_record = {}
    score_record2 = {}
    for lab in Counter(combined_labeled_array.flatten()):
        if lab == 0:
            continue  # Skip background label

        # Create a mask for the current segment
        segment_mask = (combined_labeled_array == lab).astype(np.uint8)

        # Resize the segment mask to match the score dimensions
        resized_mask = cv2.resize(segment_mask, (score.shape[1], score.shape[0]), interpolation=cv2.INTER_AREA)

        # Copy the score array and apply the mask
        each_score = score.copy()
        each_score[resized_mask < 1] = np.nan
        score_record[lab] = np.nanmean(each_score)

        # Filter out nan values and store only the non-nan scores
        non_nan_scores = each_score[~np.isnan(each_score)]
        score_record2[lab] = non_nan_scores
    # Store the score record for this file
    all_score_records[file_path] = score_record
    all_score_records2[file_path] = score_record2


# Extracting file names from paths and using them as headers
headers = ['Label'] + [os.path.splitext(os.path.basename(path))[0] for path in pickle_files]

# Preparing data for CSV
csv_data = []

# Collect all labels across files
all_labels = set()
for score_record in all_score_records.values():
    all_labels.update(score_record.keys())

# Organize data by labels
for label in sorted(all_labels):
    row = [label]
    for file_path in pickle_files:
        score_record = all_score_records[file_path]
        row.append(score_record.get(label, 'N/A'))  # 'N/A' if label not in this file
    csv_data.append(row)

# Write data to CSV file
csv_filename = f'{out_path}combined_score_records.csv'
with open(csv_filename, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(headers)
    writer.writerows(csv_data)

