"""
Handy functions to load datasets
version 1.3

Author: Matteo Pariset

"""

from unicodedata import name
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
import re

import os

from .generic import *

def as_full_dataset_name(name, embeddings_model_name, compression_strategy=None, extension='npz'):
    if compression_strategy is not None:
        name = f"{name}_{embeddings_model_name}_{compression_strategy}"
    else:
        name = f"{name}_{embeddings_model_name}"

    if extension is not None:
        return f"{name}.{extension}"
    else:
        return name

slice_pattern = re.compile('(?<=p)([0-9]+)(?=\.)')

def find_dataset_files(dataset_name, embeddings_model_name, compression_scheme, is_split=True, embs_root_name=None):
    if embs_root_name is not None:
        raise ValueError("Nested datasets not supported by this version")
    
    dir_content = os.listdir(get_dataset_dir(dataset_name))
    return filter_dataset_files(dir_content, dataset_name, embeddings_model_name, compression_scheme, is_split)

def filter_dataset_files(dir_content, dataset_name, embeddings_model_name, compression_scheme, is_split=True):
    if is_split:
        name_filter = re.compile(f'.*{dataset_name}_{embeddings_model_name}_{compression_scheme}_p[0-9]+.npz$')
    else:
        name_filter = re.compile(f'.*{dataset_name}_{embeddings_model_name}_{compression_scheme}.npz$')

    results = list(filter(name_filter.match, dir_content))
    
    # Reorder slices
    ordered_results = list(map(lambda x: int(slice_pattern.search(x)[0]), results))
    if len(ordered_results) > 0:
        reordering_idxs = np.argsort(ordered_results)
        return np.array(results)[reordering_idxs]
    else:
        return np.array([], dtype=np.string_)


def shift_and_normalize_labels(labels, ref_infos):
    # Median alignment
    shift = labels.median() - ref_infos['norm']['median']
    obt_labels = labels - shift
    # Normalization
    obt_labels -= ref_infos['norm']['mean']
    obt_labels /= ref_infos['norm']['max']

    return obt_labels

def load_centroid_mask(dataset_name, path="views/clusters/", cluster_method="usearch"):
    num_mask = pd.read_csv(relative_path(path, cluster_method, f"{dataset_name}_mask.csv"))['centroid_mask'].to_numpy()
    return num_mask

def load_cluster_mask(dataset_name, path="views/clusters/", cluster_method="usearch"):
    cluster_mask = pd.read_csv(relative_path(path, cluster_method, f"{dataset_name}_mask.csv"))['cluster_mask'].to_numpy()
    return cluster_mask

def is_centroid(centroid_mask):
    """
    Transforms centroid mask into bool array, indicating which sequences are centroids
    """
    return centroid_mask >= 0

def load_masked_dataset(embeddings_model_name, dataset_name, labels_name, compression=None, preview=True, normalize=None, centroid_mask=None, cluster_mask=None):
    embedding_files = np.char.array(os.listdir(get_dataset_dir(dataset_name)))
    filename_stub = as_full_dataset_name(dataset_name, embeddings_model_name, compression, extension=None)

    if normalize is None:
        if labels_name in ['selection_class']:
            normalize = False
        else:
            normalize = True
        print("[INFO] Automatically determining if normalization is needed. Decision:", normalize)

    # Identify desired file
    desired_file = (embedding_files.find(filename_stub) != -1) * np.logical_not(embedding_files.find('dvc') != -1)
    assert not (desired_file.sum() < 1), f"No file compatible with '{filename_stub}' found"
    
    sorted_embs_filenames = np.sort(embedding_files[desired_file])
    assert sorted_embs_filenames.shape[0] == desired_file.sum(), "[ERROR] Unable to sort slice filenames"

    if (desired_file.sum() > 1):
        print(f"[INFO] Found {desired_file.sum()} incremental savings. Reconstructing embeddings")

    unmasked_embs = np.concatenate([np.load(os.path.join(get_dataset_dir(dataset_name), f"{slice_file}"))['arr_0'] for slice_file in sorted_embs_filenames])
    unmasked_labels_df = pd.read_csv(os.path.join(get_dataset_dir(dataset_name), f"{dataset_name}.csv")).reset_index()

    if centroid_mask is None:
        embs = unmasked_embs
        labels_df = unmasked_labels_df
    else:
        print("[INFO] unclustered\t: labels_shape", unmasked_labels_df.shape)
        embs = unmasked_embs[is_centroid(centroid_mask)]
        labels_df = unmasked_labels_df.iloc[centroid_mask >= 0]  # np.repeat(np.expand_dims(mask, 1), len(labels_df.columns), axis=1))
        print("[INFO] centroids\t: labels_shape", labels_df.shape)
    
    labels = labels_df[labels_name].copy()

    infos = {
        'name': f"{dataset_name}_{labels_name}",
        'mask': centroid_mask,
        'extras': unmasked_labels_df,
    }
    
    if cluster_mask is not None:
        print("[INFO] Recording clusters: this can increase the memory footprint")
        infos.update({
            'cluster_embs': {centroid: unmasked_embs[cluster_mask == centroid] for centroid in np.unique(cluster_mask) if (cluster_mask == centroid).sum() > 1},
            'cluster_labels': {centroid: unmasked_labels_df.iloc[cluster_mask == centroid] for centroid in np.unique(cluster_mask)}
        })

    if preview:
        masked_labels = np.ma.masked_invalid(labels)
        plt.plot(np.arange(labels.shape[0])[~masked_labels.mask], masked_labels.compressed())
        plt.hlines(labels.median(), xmin=0, xmax=labels.shape[0], color='orange')
        plt.title("affinity");
        plt.show()

    if isinstance(normalize, bool) and normalize:
        orig_median = labels.median()
        orig_mean = labels.mean()
        labels -= orig_mean
        orig_max =  np.abs(labels).max()
        labels /= orig_max
        infos['norm'] = {'mean': orig_mean, 'median': orig_median, 'max': orig_max}
    elif isinstance(normalize, dict):
        labels = shift_and_normalize_labels(labels, normalize)
    else:
        infos['norm'] = {'mean': labels.mean(), 'median': labels.median(), 'max': np.abs(labels - labels.mean()).max()}

    assert embs.shape[0] == labels.shape[0], "Mismatch between dims of embeddings and labels"

    return embs.astype('float32'), labels, infos


def load_dataset(embeddings_model_name, dataset_name, labels_name, sequence_names=["sequence"], compression=None, is_split=True, preview=True, normalize=None, labels_transform=None, fillna=None, slices_cap=None):
    filename_stub = as_full_dataset_name(dataset_name, embeddings_model_name, compression, extension=None)

    if normalize is None:
        if labels_name in ['selection_class']:
            normalize = False
        else:
            normalize = True
        print("[INFO] Automatically determining if normalization is needed. Decision:", normalize)

    # Identify desired file
    sorted_embs_filenames = find_dataset_files(dataset_name, embeddings_model_name, compression_scheme=compression, is_split=is_split)
    assert len(sorted_embs_filenames) >= 1, f"No file compatible with '{filename_stub}' found"
    
    if (sorted_embs_filenames.shape[0] > 1):
        print(f"[INFO] Found {sorted_embs_filenames.shape[0]} incremental savings. Reconstructing embeddings")

    labels_df = pd.read_csv(os.path.join(get_dataset_dir(dataset_name), f"{dataset_name}.csv")).reset_index()
    if labels_name is not None:
        if fillna is not None:
            labels_df[labels_name] = labels_df[labels_name].fillna(fillna)

        labels = labels_df[labels_name].copy()

        if labels_transform is not None and callable(labels_transform):
            print(f"[INFO] Applying labels transform")
            labels = labels.apply(labels_transform)

    embs = {}
    for seq_name in sequence_names:
        print(f"[INFO] Loading sequence: '{seq_name}'")
        cursor = 0
        for slice_num, slice_file in enumerate(tqdm(sorted_embs_filenames[:slices_cap].tolist())):
            if len(sequence_names) > 1 or sequence_names[0] != "sequence":
                temp_slice = np.load(os.path.join(get_dataset_dir(dataset_name), f"{slice_file}"))[seq_name].astype('float32')
            else:
                temp_slice = np.load(os.path.join(get_dataset_dir(dataset_name), f"{slice_file}"))['arr_0'].astype('float32')

            if seq_name not in embs:
                # Instantiate entire array (to avoid doubling when copying)
                embs[seq_name] = np.zeros((temp_slice.shape[0] * sorted_embs_filenames[:slices_cap].shape[0], *temp_slice.shape[1:]), dtype="float32")

            embs[seq_name][cursor:cursor+temp_slice.shape[0]] = temp_slice[:]
            cursor += temp_slice.shape[0]
            
            del temp_slice

        embs[seq_name] = embs[seq_name][:cursor]

        if labels_name is None:
            labels = pd.Series([0] * embs[seq_name].shape[0])

        assert embs[seq_name].shape[0] == labels.shape[0], f"Mismatch between dims of embeddings and labels: {embs[seq_name].shape[0]} vs {labels.shape[0]}"

    if labels_name is not None and preview:
        masked_labels = np.ma.masked_invalid(labels)
        plt.plot(np.arange(labels.shape[0])[~masked_labels.mask], masked_labels.compressed())
        plt.hlines(labels.median(), xmin=0, xmax=labels.shape[0], color='orange')
        plt.title("affinity");
        plt.show()

    infos = {
        'name': f"{dataset_name}_{labels_name}",
        'extras': labels_df,
    }

    if isinstance(normalize, str) and normalize == "minmax":
        orig_median = labels.median()
        orig_mean = labels.mean()
        labels -= orig_mean
        orig_max =  np.abs(labels).max()
        labels /= orig_max
        infos['norm'] = {'mean': orig_mean, 'median': orig_median, 'max': orig_max}
    elif isinstance(normalize, str) and normalize == "standard":
        orig_mean = labels.mean()
        labels -= orig_mean
        orig_std = labels.std()
        labels /= orig_std
        infos['norm'] = {'mean': orig_mean, 'std': orig_std}
    elif isinstance(normalize, dict):
        labels = shift_and_normalize_labels(labels, normalize)
    else:
        infos['norm'] = {'mean': labels.mean(), 'median': labels.median(), 'max': np.abs(labels - labels.mean()).max()}


    return embs, labels, infos 


def reduce_embs(embs, strategy='mean'):
    if len(embs.shape) == 2:
        print("[WARNING] Embeddings are already reduced: not doing anything")
        return embs

    reduction_func = None
    if isinstance(strategy, str):
        if strategy == 'mean':
            reduction_func = lambda x:  x.mean(axis=1)
        elif strategy == 'cdr3_syb':
            reduction_func = lambda x: x[:,96:115,:].mean(axis=1)
        elif strategy == 'cdr2_syb':
            reduction_func = lambda x: x[:,45:56,:].mean(axis=1)
        elif strategy == 'cdr1_syb':
            reduction_func = lambda x: x[:,24:36,:].mean(axis=1)
        elif strategy == 'all_cdrs_syb':
            reduction_func = lambda x: x[:,np.s_[24:36, 45:56, 96:115],:].mean(axis=1)
        elif strategy == 'scaffold':
            reduction_func = lambda x: x[:,35:45,:].mean(axis=1)
        elif strategy == 'nanobody_vh':
            reduction_func = lambda x: x[:,:150,:].mean(axis=1)
        elif strategy == 'cdr3_ab':
            reduction_func = lambda x: x[:,80:115,:].mean(axis=1)
        elif strategy == 'bos':
            reduction_func = lambda x: x[:,0,:]
        else:
            raise ValueError("Unknown reduction strategy")
    elif isinstance(strategy, slice):
        reduction_func = lambda x: x[:,strategy,:].mean(axis=1)

    elif callable(strategy):
        reduction_func = strategy

    return reduction_func(embs)

def get_embs_dim(embeddings_model_name):
    if embeddings_model_name == "esm" or "esm-1v" in embeddings_model_name or "esm-2" in embeddings_model_name:
        embs_dim = 1280
    elif embeddings_model_name == "protein_bert":
        embs_dim = 1562
    elif embeddings_model_name == "prot_gpt2":
        embs_dim = 1280
    else:
        raise ValueError()

    return embs_dim
