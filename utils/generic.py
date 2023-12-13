"""
Handy functions to load datasets
version 1.3

Author: Matteo Pariset

"""

import os
from os.path import join

def relative_path(*rel_path, root_name=None):
    if root_name is None:
        root_name = "developability_profiling"
    
    given_path = os.path.normpath(root_name)
    if os.path.isabs(given_path):
        # Use given absolute path
        abs_root = given_path.split(os.sep)
    else:
        # Find relative root in current path
        split_curr_path = os.path.normpath(os.getcwd()).split(os.sep)
        abs_root = split_curr_path[1:split_curr_path.index(root_name)+1]
    
    return join("/", *abs_root, *rel_path)

def get_dataset_dir(dataset_name):
    if dataset_name == "developability":
        return relative_path("data", "native")
    elif dataset_name == "developability-thera":
        return relative_path("data", "mAbs")
    elif dataset_name == "developability-humanised":
        return relative_path("data", "kymouse")
    else:
        raise ValueError(f"Dataset {dataset_name} is not available")

def info(*msg, **kwargs):
    print("[INFO]", *msg, **kwargs)

def warn(*msg, **kwargs):
    print("[WARNING]", *msg, **kwargs)

def error(*msg, **kwargs):
    print("[ERROR]", *msg, **kwargs)