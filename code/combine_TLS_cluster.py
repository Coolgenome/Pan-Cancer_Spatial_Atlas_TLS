import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from utils import read_lines, save_pickle, load_pickle, save_image

prefix = sys.argv[1]
prefix2 = sys.argv[2]
prefix3 = sys.argv[3]
prefix4 = sys.argv[4]

numeric_prefix3 = int(prefix3)
numeric_prefix4 = int(prefix4)

def plot_label_masks_new(labels, prefix, names=None, white_background=True):
    cmap = plt.get_cmap('tab10')
    color_pos = 1.0
    color_neg = 0.0
    color_bg = 0.0  # black
    if white_background:
        color_bg = 0.0  # white
    foreground = labels >= 0
    labs_uniq = np.unique(labels)
    labs_uniq = labs_uniq[labs_uniq >= 0]
    if names is None:
        names = labs_uniq
    for lab in labs_uniq:
        mask = labels == lab
        img = np.zeros(mask.shape+(3,), dtype=np.float32)
        img[mask] = color_pos
        img[~mask] = color_neg
        img[~foreground] = color_bg
        img = (img * 255).astype(np.uint8)
        nam = names[lab]
        save_image(img, f'{prefix}{nam}.png')

data = load_pickle(f'{prefix}tls_cluster/nclusters00{prefix2}/labels.pickle')
data[data == numeric_prefix3] = numeric_prefix4
save_pickle(data, f'{prefix}tls_cluster/nclusters00{prefix2}/labels_combined.pickle')
plot_label_masks_new(data,f'{prefix}tls_cluster/nclusters00{prefix2}/masks_combined/')

