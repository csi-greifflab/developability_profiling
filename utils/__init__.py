"""
Some tools
version 1.4

Author: Matteo Pariset

"""

import importlib

from . import generic
from . import dataset_loader
from . import fasta
from . import graphs
from . import logs

# Useful for Pylance
from .generic import *
from .dataset_loader import *
from .fasta import *
from .graphs import *
from .logs import *

core_lib_modules = [generic, dataset_loader, fasta, graphs, logs]

def _refresh_module(module, context_module):
    module = importlib.reload(module)
    for obj in module.__dict__.keys():
        if obj[0] != "_":
            context_module.__dict__[obj] = module.__dict__[obj]

def refresh(context_module):
    # print(sys.modules)
    for m in core_lib_modules:
        _refresh_module(m, context_module)

    if importlib.util.find_spec('dask_ml') is not None:
        from . import distributed_dataset_loader

        extended_lib_modules = [distributed_dataset_loader]

        for m in extended_lib_modules:
            _refresh_module(m, context_module)

