import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from utils import load_yaml, load_pickle, write_lines, save_pickle
from visual import plot_labels, plot_colorbar
from visual import plot_matrix as plot_matrix_raw
from visual import cmap_myset
from image import upscale as upsize
from cluster import cluster, sort_labels


np.random.seed(0)


def plot_matrix(x, *args, **kwargs):
    lower, upper = 0.001, 0.999
    x = np.clip(
            x, np.nanquantile(x, lower), np.nanquantile(x, upper))
    plot_matrix_raw(
            x, *args, **kwargs)


def upscale(x, factor):
    if factor is not None and factor > 1:
        shape = np.array(x.shape) * factor
        x = upsize(x, shape)
    return x


def load_scores(structures, prefix):
    for name, content in structures.items():
        score = load_pickle(f'{prefix}{name}.pickle')
        content['score'] = score
        if 'subtypes' in content.keys():
            load_scores(content['subtypes'], prefix)


def plot_scores(scores, mask, prefix, cmap):
    scores_arr = np.array(list(scores.values()))
    minmax = np.nanmin(scores_arr), np.nanmax(scores_arr)
    score_tot = np.sum(scores_arr, axis=0)
    score_tot[~mask] = 1.0
    normalized = {}
    for name, s in scores.items():
        plot_matrix(
                s,
                f'{prefix}{name}/score-scaled.png',
                minmax=minmax, cmap=cmap,
                white_background=True)
        s = s / score_tot
        save_pickle(s, prefix+'score-scaled.pickle')
        plot_matrix(
                s,
                f'{prefix}{name}/score-normalized.png',
                minmax=(0.0, 1.0), cmap=cmap,
                white_background=True)
        normalized[name] = s
    return normalized


def draw_sample(x):
    mask = np.isfinite(x).all(-1)
    shape = x.shape[:-1]
    x = x[mask]
    r = np.random.rand(x.shape[0]) * x.sum(1)
    np.cumsum(x, axis=1)
    is_above = r[..., np.newaxis] > np.cumsum(x, axis=1)
    idx = np.sum(is_above, axis=1)
    label = np.full(shape, -1)
    label[mask] = idx
    return label


def remap_labels(labels, label_map):
    labels_new = labels.copy()
    for i in range(labels.max()+1):
        labels_new[labels == i] = label_map[i]
    return labels_new


def get_centers(labels, x):
    n_labels = labels.max() + 1
    n_features = x.shape[-1]
    centers = np.full((n_labels, n_features), np.nan)
    for i in range(n_labels):
        isin = labels == i
        centers[i] = x[isin].mean(0)
    return centers


def detect_dominance(signi, probs, threshold):
    assert threshold >= 0.5
    for i in range(probs.shape[0]):
        if (probs[i] > threshold).any():
            j = probs[i].argmax()
            signi[i] = False
            signi[i, j] = True
    return signi


def multiclass(
        x, mask, n_clusters=10, adjust_clusters=False,
        threshold_dominance=0.6, threshold_significance=0.1):

    x = np.stack(x)
    x[:, ~mask] = np.nan

    labels, __ = cluster(x, n_clusters=n_clusters, method='km')
    probs = get_centers(labels, x.transpose(1, 2, 0))

    if adjust_clusters:
        signi = probs > threshold_significance
        signi = detect_dominance(signi, probs, threshold_dominance)
        __, label_map = np.unique(
                signi, axis=0, return_inverse=True)
        labels = remap_labels(labels, label_map)
        labels, __ = sort_labels(labels)
        probs = get_centers(labels, x.transpose(1, 2, 0))

    return labels, probs


def classify(
        x, mask, prefix,
        use_multiclass=False, n_clusters=10, use_sampling=False,
        adjust_clusters=False):

    if use_multiclass:
        labels, probs = multiclass(
                x, mask, n_clusters=n_clusters,
                adjust_clusters=adjust_clusters)
        out = labels, probs
    else:
        if use_sampling:
            labels = draw_sample(np.stack(x, -1))
        else:
            labels = np.argmax(x, axis=0)
        labels[~mask] = -1
        out = labels
    return out


def plot_probs(x, filename, cat_names=None, cmap=None):
    x = np.round(x * 100).astype(int)
    df = pd.DataFrame(x)
    if cat_names is not None:
        df.columns = cat_names
    sns.heatmap(
            df, cmap='magma',
            annot=True, fmt='02d', annot_kws={'fontsize': 20},
            square=True, linewidth=0.5,
            cbar_kws=dict(label='Percentage'))
    plt.yticks(rotation=0)
    plt.tick_params(
            axis='both', which='major',
            bottom=False, labelbottom=False,
            top=False, labeltop=True,
            left=False, labelleft=True,
            right=False, labelright=False,
            )

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(filename)


def plot_cmap(
        cmap, n, filename, fmt='01d', label_inside=False, label_left=True):
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    colors = [cmap(i) for i in range(n)]
    cmap = LinearSegmentedColormap.from_list(
            'mycmap', colors, len(colors))
    df = pd.DataFrame(np.arange(n)[: np.newaxis])
    df.columns = ['Cluster']
    sns.heatmap(
            df, cmap=cmap,
            annot=label_inside, fmt=fmt, annot_kws={'fontsize': 20},
            square=True, linewidth=0.5, cbar=False)
    plt.tick_params(
            axis='both', which='major',
            bottom=False, labelbottom=False,
            top=False, labeltop=True,
            left=False, labelleft=label_left,
            right=False, labelright=False)
    if label_left:
        plt.yticks(rotation=0)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(filename)


def classify_multiclass(
        normalized, mask, prefix, cat_names,
        adjust_clusters, cmap='tab10', n_clusters=10):

    n_clusters = 10
    labels_mc, probs_mc = classify(
            normalized, mask, prefix,
            use_multiclass=True, n_clusters=n_clusters,
            adjust_clusters=adjust_clusters)
    plot_labels(
            labels_mc, cmap=cmap,
            white_background=True,
            filename=prefix+'labels-multicls.png')
    plot_probs(
            probs_mc, cat_names=cat_names, cmap=cmap,
            filename=prefix+'labels-multicls-distri.png')
    n_labels = labels_mc.max() + 1
    if n_labels > 1:
        plot_cmap(
                cmap, n=n_labels,
                filename=prefix+'labels-multicls-cmap.png')


def compare(scores, mask, prefix, cmap):
    # Existing code to generate normalized
    normalized = plot_scores(scores, mask, prefix, cmap=cmap)
    normalized = list(normalized.values())
    names = np.array(list(scores.keys()))

    # Existing code to generate cmap
    cmap = cmap_myset
    n_colors = len(names)
    plot_cmap(cmap=cmap, n=n_colors, filename=prefix+'labels-cmap.png')

    labels_maxprob = classify(normalized, mask, prefix, use_sampling=False)
    save_pickle(labels_maxprob, prefix+'labels.pickle')
    # Existing code to plot and write labels
    plot_labels(labels_maxprob, cmap=cmap, white_background=True, filename=prefix+'labels.png')
    write_lines(names, prefix+'labels-names.txt')

    # New code to generate separate figures for each label
    unique_labels = np.unique(labels_maxprob[labels_maxprob >= 0])  # Ignore labels with -1
    for label in unique_labels:
        label_mask = (labels_maxprob == label)
        single_label_image = np.full(labels_maxprob.shape, -1)  # Initialize with -1
        single_label_image[label_mask] = label
        plot_labels(single_label_image, cmap=cmap, white_background=True, filename=f'{prefix}label_{label}.png')
    


def analyze(
        structures, prefix, mask=None,
        fill_with_min=True, upscale_factor=None):
    cmap = 'turbo'
    plot_colorbar(
            cmap, n_labels=None, filename=prefix+'colorbar.png')
    score_dict = {}
    for name, content in structures.items():
        score = content['score'].copy()
        if upscale_factor is not None:
            score = upscale(score, factor=upscale_factor)
        filler_val = np.nan
        if fill_with_min:
            filler_val = np.nanmin(score)
        is_fin = np.isfinite(score)

        if mask is not None:
            score[(~mask)*is_fin] = filler_val
        plot_matrix(
                score, f'{prefix}{name}/score.png', cmap=cmap,
                white_background=True)

        score_dict[name] = score

        if 'threshold' in content.keys():
            threshold = content['threshold']
            mask_sub = score >= threshold
            if mask is not None:
                mask_sub = mask_sub * mask

            score_thresholded = score.copy()
            score_thresholded[(~mask_sub)*is_fin] = filler_val
            plot_matrix(
                    score_thresholded,
                    f'{prefix}{name}/score-thresholded.png',
                    cmap=cmap, white_background=True)
        else:
            mask_sub = mask

        if 'subtypes' in content.keys():
            score_sub_dict = analyze(
                    content['subtypes'], f'{prefix}{name}/',
                    mask=mask_sub)
            compare(
                    score_sub_dict, mask=mask_sub,
                    prefix=f'{prefix}{name}/', cmap=cmap)
    return score_dict


def main():

    prefix = sys.argv[1]  # e.g. 'data/G123/markers/phenotype/'

    structures = load_yaml('/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_signature_final/structures.yml')
    load_scores(structures, prefix+'markers/phenotype/raw/')

    analyze(structures, prefix+'analysis/')

if __name__ == '__main__':
    main()
