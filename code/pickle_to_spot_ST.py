import os
import sys
import pickle
import PIL
from PIL import Image
import numpy as np
import pandas as pd
import cv2
import csv

prefix = sys.argv[1]
pickle_file = sys.argv[2]
out_suffix = sys.argv[3]
csv_path = sys.argv[4]

PIL.Image.MAX_IMAGE_PIXELS = None

def load_image(filename, verbose=True):
    img = Image.open(filename)
    img = np.array(img)
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]  # remove alpha channel
    if verbose:
        print(f'Image loaded from {filename}')
    return img


def load_pickle(filename, verbose=True):
    with open(filename, 'rb') as file:
        x = pickle.load(file)
    if verbose:
        print(f'Pickle loaded from {filename}')
    return x

def save_as_csv(data, filename):
    data.to_csv(filename, index=True)
    print(f'Saved CSV file: {filename}')

pickle_or_fig = prefix+pickle_file

df = pd.read_csv(csv_path)

tls_ids = df[df['Manual'] == 'TLS']['TLS_ID'].tolist()

cnts = load_pickle(pickle_or_fig)
raw_img = load_image(prefix+'he.jpg')
size_ratio = raw_img.shape[0]/cnts.shape[0]
crop_dict = {}

with open(prefix+'radius-raw.txt', 'r') as file:
    radius_raw = file.read()
    radius_numeric = float(radius_raw)
    radius_numeric = radius_numeric * 2 

radius = round(radius_numeric/(2*size_ratio))
image = np.array(cnts)
locs = pd.read_csv(prefix+'locs.tsv',index_col=0,sep='\t')
locs = (locs/size_ratio).astype(np.int64)

for each,cell_index in zip(np.array(locs),locs.index):
    #Defines a square region centered at each with a side length of 2 * radius.
    y1, x1, y2, x2 = each[0]-radius, each[1]-radius, each[0]+radius, each[1]+radius
    cropped_region = image[x1:x2, y1:y2]
    #Ensures the cropped region is a square of size 2 * radius.
    if (cropped_region.shape[0] != cropped_region.shape[1]) | (cropped_region.shape[0] != radius*2):
        continue
    
    #Creates a mask of the same shape as cropped_region with an additional dimension for color channels.
    mask = np.zeros(list(cropped_region.shape)+[3])
    #Draws a filled circle in the mask with the center at (radius, radius) and a radius of radius.
    mask = cv2.circle(mask, (radius,radius), radius, (255, 255, 255), thickness=-1)
    #Converts the mask to a binary mask where pixels inside the circle are set to False.
    mask = mask.mean(2)<10
    #Sets the values in cropped_region outside the circle to zero.
    cropped_region[mask] = 0
    unique_values = np.unique(cropped_region)
    # Check unique values and set crop_dict accordingly
    if len(unique_values) == 1 and unique_values[0] == 0:
        crop_dict[cell_index] = 0
    else:
        # Exclude zero and take the first non-zero unique value
        non_zero_values = unique_values[unique_values != 0]
        if len(non_zero_values) > 0:
            crop_dict[cell_index] = non_zero_values[0]
        else:
            crop_dict[cell_index] = 0
    
spot = pd.DataFrame.from_dict(crop_dict, orient='index', columns=['TLS_Location'])
spot['TLS_Location'] = spot['TLS_Location'].apply(lambda x: x if x in tls_ids else 0)

combined_labeled_arrays_csv_path = prefix+out_suffix+'combined_labeled_array.csv'
save_as_csv(spot, combined_labeled_arrays_csv_path)
