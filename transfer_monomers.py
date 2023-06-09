'''For transferring residue-averaged charges from reductions to full-sized collection'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path

# Polymer Imports
from polysaccharide.polymer.management import PolymerManager, MolFilterBuffer
from polysaccharide.polymer.filtering import has_monomers_chgd

# Utility function imports
from workflow_functs import retrieve_monomers

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__ # use script docstring as help description 
)
parser.add_argument('-src' , '--source_name', help='The name of the collection of reduced Polymers to draw monomers from', required=True)
parser.add_argument('-targ', '--target_name', help='The name of the target output collection of Polymers to move charged monomers to', required=True)
MolFilterBuffer.argparse_inject(parser)

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

# defining paths
source_path = COLL_PATH / args.source_name
output_path = COLL_PATH / args.output_name

# defining filters
molbuf = MolFilterBuffer.from_argparse(args)
molbuf.charges = True # force preference for charges 
mol_filters = molbuf.filters
mol_filters.append(has_monomers_chgd) # also assert that, not only do charges exist, but that they've been monomer-averaged

# Execution
# ------------------------------------------------------------------------------

from logging import Logger
if __name__ == '__main__':
    # load small collection and clone into new folder (over OLD collection)
    mgr_small = PolymerManager(source_path)
    mgr_large = PolymerManager(output_path) # NOTE : load after cloning to ensure collection is updated

    transfer_chgd_mono = mgr_large.logging_wrapper(
        proc_name='Charged monomer transfer',
        filters=mol_filters
    )(retrieve_monomers)

    transfer_chgd_mono(mgr_small)