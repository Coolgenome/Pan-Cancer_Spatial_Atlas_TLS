import sys
import PIL
from PIL import Image, ImageDraw, ImageFont, ImageColor
import numpy as np
#from scipy.ndimage import label 
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

# Concatenate to create the full paths
mask_path= prefix + mask_suffix
out_path = prefix + out_suffix
os.makedirs(out_path, exist_ok=True)

PIL.Image.MAX_IMAGE_PIXELS = None

# Load your image and mask
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

score_he = Image.open(prefix + 'analysis/lymphoid/score_resized.png')
cd3d_img = Image.open(prefix + 'cnts-super-plots/CD3D.png')
cd20_img = Image.open(prefix + 'cnts-super-plots/MS4A1.png')
cr2_img = Image.open(prefix + 'cnts-super-plots/CR2.png')
fcer2_img = Image.open(prefix + 'cnts-super-plots/FCER2.png')
mki67_img = Image.open(prefix + 'cnts-super-plots/MKI67.png')
bcl6_img = Image.open(prefix + 'cnts-super-plots/BCL6.png')

cd3d_img_resized = cd3d_img.resize(score_he.size)
cd20_img_resized = cd20_img.resize(score_he.size)
cr2_img_resized = cr2_img.resize(score_he.size)
fcer2_img_resized = fcer2_img.resize(score_he.size)
mki67_img_resized = mki67_img.resize(score_he.size)
bcl6_img_resized = bcl6_img.resize(score_he.size)

df = pd.read_csv(csv_path)
tls_ids = set(df['TLS_ID'])

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


def process_image(img, df, output_filename):
    # Define color maps
    color_map = {
        'TLS': '#0a3b03',  # Greenish color
        'TO-P': '#433578',  # Purple color
        # Add more mappings as needed for other possible values
    }
    color_map2 = {
    'TLS': '#ede9f2',  # Purple white
    'TO-P': '#f7b5a6',  # Salmon
    }  

    manual_map = dict(zip(df['TLS_ID'], df['Manual']))
    filter_map = dict(zip(df['TLS_ID'], df['Filter']))
    # Select color map based on output_filename
    selected_color_map = color_map2 if output_filename.startswith('2_') or output_filename.startswith('5_') else color_map
    offset = 5
    # Convert the original PIL image to a NumPy array for contour drawing
    img_plot_np = np.array(img.convert('RGB'))

    for i in range(1, combined_num_features + 1):
        if i in tls_ids:
            manual_class = manual_map.get(i)
            color = selected_color_map.get(manual_class, 'black')

            segment_mask = (combined_labeled_array == i).astype(np.uint8) * 255
            segment_mask_cv = np.array(Image.fromarray(segment_mask))

            contours, _ = cv2.findContours(segment_mask_cv, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

            if contours:
                color_rgb = ImageColor.getrgb(color)
                color_bgr = color_rgb[::-1]
                # Draw contours on the NumPy array
                cv2.drawContours(img_plot_np, contours, -1, color_bgr, 20)

    # Convert the NumPy array with drawn contours back to a PIL image
    img_with_contours = Image.fromarray(img_plot_np)
    draw = ImageDraw.Draw(img_with_contours)
    font_size = 225
    font = ImageFont.truetype("/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/Arial.ttf", font_size)

    used_positions = []
    for i in range(1, combined_num_features + 1):
        if i in tls_ids:
            manual_class = manual_map.get(i)
            color = selected_color_map.get(manual_class, 'black')

            segment_mask = (combined_labeled_array == i).astype(np.uint8) * 255
            segment_mask_cv = np.array(Image.fromarray(segment_mask))
            contours, _ = cv2.findContours(segment_mask_cv, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

            if contours:
                x, y, w, h = cv2.boundingRect(contours[0])
                text_position = (x, max(y - font_size - offset, 0))
                text_to_draw = f"{i}"
                #text_width, text_height = font.getsize(text_to_draw)
                bbox = font.getbbox(text_to_draw)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]

                for used_position in used_positions:
                    if abs(text_position[1] - used_position[1]) < text_height + offset:
                        text_position = (text_position[0], used_position[1] + text_height + offset)
                        break

                used_positions.append(text_position)

                # Draw the index number on the image
                draw.text(text_position, text_to_draw, fill=color, font=font)
            else:
                print(f"Segment {i} not found or is empty.")

    # Save the final image
    img_with_contours.save(f'{out_path}{output_filename}')


def process_bw_image(df, output_filename='3_BW_highlighted_image.png'):

    color_map = {
        'TLS': '#FFFFFF',
        'TO-P': '#433578',
    }

    manual_map = dict(zip(df['TLS_ID'], df['Manual']))
    image_width, image_height = he_image.size

    font_size = 225
    font = ImageFont.truetype("/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/Arial.ttf", font_size)

    used_positions = []
    black_background = np.zeros((image_height, image_width), dtype=np.uint8)
    text_regions_mask = np.zeros((image_height, image_width), dtype=np.uint8)

    # First pass: collect positions and white mask regions
    for i in range(1, combined_num_features + 1):
        if i in tls_ids:
            segment_mask = (combined_labeled_array == i).astype(np.uint8) * 255
            bbox = Image.fromarray(segment_mask).getbbox()
            if bbox:
                text_regions_mask[combined_labeled_array == i] = 255

                offset = 5
                text_position = (bbox[0], max(bbox[1] - font_size - offset, 0))
                text_to_draw = f"{i}"
                bbox = font.getbbox(text_to_draw)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]

                for used_position in used_positions:
                    if abs(text_position[1] - used_position[1]) < text_height + offset:
                        text_position = (text_position[0], used_position[1] + text_height + offset)
                        break

                if text_position[0] + text_width > image_width:
                    text_position = (image_width - text_width, text_position[1])
                if text_position[1] < 0:
                    text_position = (text_position[0], 0)

                used_positions.append(text_position)

    result_image = np.where(text_regions_mask[..., np.newaxis], 255, black_background[..., np.newaxis])
    result_image_2d = result_image.reshape((image_height, image_width))
    result_image_pil = Image.fromarray(result_image_2d)
    draw = ImageDraw.Draw(result_image_pil)

    # Second pass: draw index text
    for i in range(1, combined_num_features + 1):
        if i in tls_ids:
            manual_class = manual_map.get(i)
            color = color_map.get(manual_class, 'black')

            segment_mask = (combined_labeled_array == i).astype(np.uint8) * 255
            bbox = Image.fromarray(segment_mask).getbbox()
            if bbox and used_positions:
                draw.text(used_positions.pop(0), f"{i}", fill=color, font=font)

    result_image_pil.save(f'{out_path}{output_filename}')


process_image(he_image,df,'1_HE_highlighted_image.png')
process_image(score_he,df,'2_TLS_score_highlighted_image.png')
process_image(cd3d_img_resized,df,'4_CD3D_highlighted_image.png')
process_image(cd20_img_resized,df,'5_MS4A1_highlighted_image.png')
process_image(cr2_img_resized,df,'6_CR2_highlighted_image.png')
process_image(fcer2_img_resized,df,'7_FCER2_highlighted_image.png')
process_image(mki67_img_resized,df,'8_MKI67_highlighted_image.png')
process_image(bcl6_img_resized,df,'9_BCL6_highlighted_image.png')
process_bw_image(df,'3_BW_highlighted_image.png')