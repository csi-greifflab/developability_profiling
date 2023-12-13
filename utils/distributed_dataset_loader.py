"""
Handles distributed datasets
version 1.0

Author: Matteo Pariset

"""

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import ast
from textwrap import wrap

import matplotlib as mpl
import matplotlib.pyplot as plt

import sklearn as sk
from dask_ml.decomposition import PCA
from umap.umap_ import UMAP

import tensorflow as tf
from tensorboard.plugins import projector

import os
from os.path import join

import multiprocessing as mp

import traceback


from dask.distributed import LocalCluster, Client
from dask_cloudprovider.aws import FargateCluster

from s3fs import S3FileSystem

import dask as dk
import dask.array as da

from configparser import ConfigParser


from .dataset_loader import *

# TODO: Handle variable slice sizes!
first_slice_size = 2589
std_slice_size = 2560
last_slice_size = 1560
embs_seq_len = 175

def load_distributed_dataset(
        dataset_name,
        compression,
        preview=False,
        normalize=None,
        centroid_mask=None,
        cluster_mask=None,
        cfg_filepath=None,
        embeddings_model_name=None
    ):

    # Load config
    cfg_parser = ConfigParser()
    cfg_parser.read(cfg_filepath);

    aws_key = cfg_parser['auth']['aws_key']
    aws_secret = cfg_parser['auth']['aws_secret']

    @dk.delayed
    def lazy_npz_load(file_path):
        from s3fs.core import S3FileSystem
        s3 = S3FileSystem(key=aws_key, secret=aws_secret)
        return np.load(s3.open(file_path))['arr_0']

    s3 = S3FileSystem(key=aws_key, secret=aws_secret, default_fill_cache=False)

    embedding_files = np.char.array(s3.listdir("adaptyv-ml-binding/" + dataset_name + "/embeddings", detail=False)[1:])
    
    if normalize is not None:
        raise NotImplementedError("Automatic label processing is not yet available on distributed datasets")

    # Identify desired file
    sorted_embs_filenames = filter_dataset_files(embedding_files, dataset_name, embeddings_model_name, compression, is_split=True)
    assert sorted_embs_filenames.shape[0] >= 1, "[ERROR] No compatible embedding files found"

    info(f"Found {sorted_embs_filenames.shape[0]} compatible slices")

    # unmasked_embs_0 = da.from_array(np.load(s3.open(sorted_embs_filenames[0]))['arr_0'])
    # unmasked_embs_reg_sl = da.concatenate([unmasked_embs_0] + [da.from_delayed(lazy_npz_load(f_name), shape=(slice_size, unmasked_embs_0.shape[1], get_embs_dim(embeddings_model_name)), dtype='float32') for f_name in sorted_embs_filenames[1:-1]])
    # last_slice = da.from_array(np.load(s3.open(sorted_embs_filenames[-1]))['arr_0'])
    # unmasked_embs = da.concatenate([unmasked_embs_reg_sl, last_slice])

    unmasked_embs_0 = da.from_delayed(lazy_npz_load(sorted_embs_filenames[0]), shape=(first_slice_size, embs_seq_len, get_embs_dim(embeddings_model_name)), dtype='float32') 
    unmasked_embs_reg_sl = da.concatenate([unmasked_embs_0] + [
        da.from_delayed(lazy_npz_load(f_name), shape=(std_slice_size, embs_seq_len, get_embs_dim(embeddings_model_name)), dtype='float32') for f_name in sorted_embs_filenames[1:-1]
    ])
    last_slice = da.from_delayed(lazy_npz_load(sorted_embs_filenames[-1]), shape=(last_slice_size, embs_seq_len, get_embs_dim(embeddings_model_name)), dtype='float32') 
    unmasked_embs = da.concatenate([unmasked_embs_reg_sl, last_slice])

    info("Distributed array successfully initialized")
    
    if centroid_mask is not None or cluster_mask is not None:
        raise NotImplementedError("Clustering is not yet supported on distributed datasets")

    return unmasked_embs


def launch_cluster(cfg_filepath, workers_num = 45):
    # Load config
    cfg_parser = ConfigParser()
    cfg_parser.read(cfg_filepath);

    info("Creating cluster...")
    cluster = FargateCluster(
        fargate_spot=True,
        worker_cpu=4096,
        worker_mem=30720,
        cluster_name_template='developability-dask-{uuid}',
        region_name=cfg_parser['cluster']['region'],
        # vpc=cfg_parser['cluster']['vpc'],
        scheduler_timeout='60 minutes',
        # security_groups=[cfg_parser['cluster']['base_security_group']],
        environment={'EXTRA_PIP_PACKAGES': 's3fs'},
        n_workers=1
    )
    info("Cluster infrastructure ready")

    client = Client(cluster)
    info("Dashboard (watch computation unfold):", client.dashboard_link)

    info(f"Spawning {workers_num} workers...")
    cluster.scale(workers_num)
    client.wait_for_workers(workers_num)
    info("Workers ready")

    return cluster, client
