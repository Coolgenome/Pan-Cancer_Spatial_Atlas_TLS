from PIL import Image
import sys
import os

prefix = sys.argv[1]

def resize_image(input_path, output_path, scale_factor):
    # Open the image
    img = Image.open(prefix+input_path)

    # Calculate the new dimensions
    new_width = int(img.width * scale_factor)
    new_height = int(img.height * scale_factor)

    # Resize the image
    resized_img = img.resize((new_width, new_height), Image.LANCZOS)

    # Save the resized image
    resized_img.save(prefix+output_path)

# Define your parameters
input_path1 = 'analysis/lymphoid/score.png'  # Path to your downscaled image
output_path1 = 'analysis/lymphoid/score_resized.png'  # Path to save the upscaled image
scale_factor = 16  # The factor by which you want to scale up

# Resize the image
resize_image(input_path1, output_path1, scale_factor)