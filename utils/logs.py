"""
Logging utils
version 1.0

Author: Matteo Pariset

"""

import numpy as np
import pandas as pd
import os

import ast

def read_last_line(filepath):
    with open(filepath, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        last_line = f.readline().decode()
    return last_line

def update_log_file(filepath, idx, *query):
    with open(filepath, 'a') as f:
        f.write(f"{idx},{','.join(query)}\n")

def init_file(file_dir, file_name, init_content):
    if not os.path.exists(file_dir):
        os.mkdir(file_dir)

    log_filepath = os.path.join(file_dir, file_name)
    if not os.path.exists(log_filepath):
        with open(log_filepath, 'w') as f:
            f.write(init_content)

    return log_filepath

def log_info(info_file, query_idx, content, logfile_dir="logs/"):
    log_filepath = init_file(logfile_dir, info_file, "")
    update_log_file(log_filepath, query_idx, content)


