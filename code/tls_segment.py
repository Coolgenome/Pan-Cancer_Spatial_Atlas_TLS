import sys
import PIL
from PIL import Image, ImageDraw, ImageFont
import numpy as np
#from scipy.ndimage import label 
import scipy.ndimage as ndi
import pickle
from collections import Counter
import cv2
import csv
import os
from utils import load_image, save_image

prefix = sys.argv[1]
mask_suffix = sys.argv[2]
out_suffix = sys.argv[3]

sig_suffix = "markers/phenotype/raw/"
gene_suffix = "cnts-super/"

# Concatenate to create the full paths
mask_path= prefix + mask_suffix
sig_pickle_path = prefix + sig_suffix
gene_pickle_path = prefix + gene_suffix
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

###image with TLS regions
filename_whole = f'{out_path}image_with_segmented_regions.png'
save_image(he_image_np,filename_whole)

mask_tls_array = np.array(mask_tls_resized)
labeled_array, num_features = ndi.label(mask_tls_array > 0)

# Loop throhe_image each feature/segment
segmented_images = []
for i in range(1, num_features + 1):
    # Create a mask for the current segment
    segment_mask = (labeled_array == i).astype(np.uint8) * 255

    # Convert the NumPy array back to PIL image for cropping
    segment_mask_pil = Image.fromarray(segment_mask)

    # Apply this mask to the image
    he_image.putalpha(segment_mask_pil)

    # Get the bounding box and crop the image
    bbox = segment_mask_pil.getbbox()
    cropped_segment = he_image.crop(bbox)

    # Store or process the cropped segment
    segmented_images.append(cropped_segment)

    # Save each segment with a unique filename
    filename = f'{out_path}segmented_image_{i}.png'
    cropped_segment.save(filename)

# Create a drawing context
he_image_cropped = Image.open(f'{out_path}image_with_segmented_regions.png')
draw = ImageDraw.Draw(he_image_cropped)
font_size = 255
font = ImageFont.truetype("/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/Arial.ttf", font_size)

# Loop throhe_image each feature/segment to draw indexes
for i in range(1, num_features + 1):
    # Create a mask for the current segment
    segment_mask = (labeled_array == i).astype(np.uint8) * 255

    # Convert the NumPy array back to PIL image for cropping
    segment_mask_pil = Image.fromarray(segment_mask)

    # Get the bounding box of the current segment
    bbox = segment_mask_pil.getbbox()

    # Determine a position for the index text (e.g., top-left corner of the bbox)
    text_position = (bbox[0], bbox[1])

    # Draw the index number on the image
    draw.text(text_position, str(i), fill='red',font=font)  # Change 'red' to your preferred text color

# Save the image with indexes
he_image_cropped.save(f'{out_path}image_with_segmented_regions_and_indexes.png')

###signature
# List of pickle file paths

pickle_files = [
    f'{sig_pickle_path}tls.pickle',
]
# Dictionary to store the loaded data from each file
loaded_data = {}

# Iterate over the list of file paths and load each file
for file_path in pickle_files:
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
        loaded_data[file_path] = data

# Assuming labeled_array is the array with labeled regions from previous steps
# and mask_tls_array is the mask array

# Dictionary to store score records for each file
all_score_records = {}

# Process each set of scores
for file_path, score in loaded_data.items():
    score_record = {}
    for lab in Counter(labeled_array.flatten()):
        if lab == 0:
            continue  # Skip background label

        # Create a mask for the current segment
        segment_mask = (labeled_array == lab).astype(np.uint8)

        # Resize the segment mask to match the score dimensions
        resized_mask = cv2.resize(segment_mask, (score.shape[1], score.shape[0]), interpolation=cv2.INTER_AREA)

        # Copy the score array and apply the mask
        each_score = score.copy()
        each_score[resized_mask < 1] = np.nan
        score_record[lab] = np.nanmean(each_score)

    # Store the score record for this file
    all_score_records[file_path] = score_record

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

###image with TLS score
# Load your image and mask
score_he = Image.open(prefix + 'analysis/lymphoid/score_resized.png')
mask_tls_for_score = Image.open(mask_path).convert('L')  # Convert mask to grayscale
mask_tls_for_score_resized = mask_tls_for_score.resize(score_he.size)

score_he_np = np.array(score_he)
mask_tls_for_score_np = np.array(mask_tls_for_score_resized)
    
mask_tls_for_score_np = mask_tls_for_score_np > 0
score_he_np[~mask_tls_for_score_np] = 0

filename_whole = f'{out_path}score_with_segmented_regions.png'
save_image(score_he_np,filename_whole)

# Create a drawing context
draw = ImageDraw.Draw(score_he)
font = ImageFont.truetype("/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/Arial.ttf", 20)

# Label connected components
mask_tls_array = np.array(mask_tls_for_score_resized)
labeled_array, num_features = ndi.label(mask_tls_array > 0)

# Process each segment
segmented_images = []
for i in range(1, num_features + 1):
    # Create a mask for the current segment
    segment_mask = (labeled_array == i).astype(np.uint8) * 255

    # Convert the NumPy array back to PIL image for cropping
    segment_mask_pil = Image.fromarray(segment_mask)

    # Apply this mask to the image
    score_he.putalpha(segment_mask_pil)

    # Get the bounding box of the current segment
    bbox = segment_mask_pil.getbbox()
    cropped_segment = score_he.crop(bbox)

    # Store or process the cropped segment
    segmented_images.append(cropped_segment)

    # Save each segment with a unique filename
    filename = f'{out_path}score_segmented_image_{i}.png'
    cropped_segment.save(filename)

# Create a drawing context
score_cropped = Image.open(f'{out_path}score_with_segmented_regions.png')
draw = ImageDraw.Draw(score_cropped)
font_size = 255
font = ImageFont.truetype("/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/Arial.ttf", font_size)

# Loop throhe_image each feature/segment to draw indexes
for i in range(1, num_features + 1):
    # Create a mask for the current segment
    segment_mask = (labeled_array == i).astype(np.uint8) * 255

    # Convert the NumPy array back to PIL image for cropping
    segment_mask_pil = Image.fromarray(segment_mask)

    # Get the bounding box of the current segment
    bbox = segment_mask_pil.getbbox()

    # Determine a position for the index text (e.g., top-left corner of the bbox)
    text_position = (bbox[0], bbox[1])

    # Draw the index number on the image
    draw.text(text_position, str(i), fill='red',font=font)  # Change 'red' to your preferred text color

# Save the image with indexes
score_cropped.save(f'{out_path}score_with_segmented_regions_and_indexes.png')


###number of pixels
segment_pixel_counts = {}

for i in range(1, num_features + 1):
    # Create a mask for the current segment
    segment_mask = (labeled_array == i).astype(np.uint8) * 255

    # Calculate the number of pixels in the current segment
    num_pixels = np.count_nonzero(segment_mask)

    # Store the pixel count in the dictionary
    segment_pixel_counts[i] = num_pixels

csv_pixel_filename = f'{out_path}segment_pixel_counts.csv'

with open(csv_pixel_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Writing the header
    writer.writerow(['Label', 'Pixel Count'])

    # Writing the data
    for segment_label, pixel_count in segment_pixel_counts.items():
        writer.writerow([segment_label, pixel_count])

