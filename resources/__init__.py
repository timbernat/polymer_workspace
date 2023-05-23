from pathlib import Path
from typing import Iterable

import pkgutil, importlib
import importlib.resources as impres
import logging
LOGGER = logging.getLogger(__name__)

def non_dunder(path : Path) -> Iterable[str]:
    '''Return all subpaths of a directory path which contain no double underscores'''
    assert(path.is_dir())
    return [file.name 
        for file in path.iterdir()
            if '__' not in file.name
    ]

AVAIL_RESOURCES = {} # load submodules, record available file assets
for _loader, _module_name, _ispkg in pkgutil.iter_modules(__path__):
    module = importlib.import_module(f'{__package__}.{_module_name}')
    globals()[_module_name] = module # register module to namespace
    AVAIL_RESOURCES[_module_name] = non_dunder(impres.files(module))

