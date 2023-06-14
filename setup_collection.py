'''(Re)populate a collection of Polymer directories from a source directory of structural and monomer information'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler

# Generic imports
import argparse
from pathlib import Path

# Custom Imports
from polysaccharide.polymer.management import PolymerManager

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__
)

parser.add_argument('-pdb', '--struct_input', help='The name of the directory to source PDB structure from', type=Path, required=True)
parser.add_argument('-mono' , '--mono_input', help='The name of the directory to source JSON monomer files from', type=Path)
parser.add_argument('-out', '--output_path' , help='The name of the output collection of Polymers', type=Path, required=True)

parser.add_argument('-r'  , '--reset'       , help='If set, will delete the target collection if it already exists'        , action='store_true')
parser.add_argument('-ps' , '--purge_sims'  , help='If set, will delete any extant MD simulations in the target collection', action='store_true')
parser.add_argument('-pl' , '--purge_logs'  , help='If set, will delete any extant log files in the target collection'     , action='store_true')

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

## defining paths
structure_path   = args.struct_input
monomer_path     = args.mono_input
collection_path  = args.output_path

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(collection_path)

    # Purge logs (if desired), then set up new loggesr
    if args.purge_logs: # NOTE : must be done BEFORE log FileHandler is created, as this will destroy it's output as well
        mgr.purge_logs(really=True)

    # Perform manager file setup / purge actions
    with ProcessLogHandler(filedir=mgr.log_dir, loggers=LOGGERS_MASTER, proc_name=f'Setup of {mgr.collection_dir.name}', timestamp=True):
        if args.reset:
            mgr.purge_collection(really=True, purge_logs=False) # Explicitly DON'T purge logs here (will be done prior to entering log loop)

        if args.purge_sims:
            mgr.purge_sims(really=True)

        if not mgr.polymers: # will be empty if not yet instantiated or if reset prior
            mgr.populate_collection(struct_dir=structure_path, monomer_dir=monomer_path)