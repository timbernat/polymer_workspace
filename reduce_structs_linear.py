'''Creates directory of reduced-chain-length structure PDBs and monomer JSONs from a collection of linear Polymers'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path

# Polymer Imports
from polysaccharide.polymer.management import PolymerManager
from polysaccharide.polymer.filtering import is_unsolvated

# Utility function imports
from workflow_functs import generate_reduced_pdb

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__ # use script docstring as help description 
)
parser.add_argument('-src', '--source_name'    , help='The name of the target collection of Polymers', required=True)
parser.add_argument('-out', '--struct_output'  , help='The name of the directory to output generated structures to (within the set PDB folder)')
parser.add_argument('-N', '--max_chain_len'    , help='Maximum number of atoms in any of the reduced chain generated. If this is specified, CANNOT specify DOP', type=int)
parser.add_argument('-D', '--DOP'              , help='The number of monomer units to include in the generated reductions.  If this is specified, CANNOT specify max_chain_len', type=int)
parser.add_argument('-lim', '--chain_len_limit', help='The maximum allowable size for a chain to be built to; any chains attempted to be built larger than this limit will raise an error', type=int, default=300)
parser.add_argument('-f', '--flip_term_labels' , help='Names of the chains on which to reverse the order of head/tail terminal group labels (only works for linear homopolymers!)', action='store', nargs='+', default=tuple())

args = parser.parse_args()
if args.struct_output is None:
    args.struct_output = f'{args.source_name}_reduced' # can't just set as generic default, since this value depends on the source input

# Arg processing
# ------------------------------------------------------------------------------

# defining paths
source_path = COLL_PATH / args.source_name

reduced_dir = COMPAT_PDB_PATH / args.struct_output # CRITICAL : Note that this is a PDB STRUCTURES path, NOT a COLLECTION path!!
reduced_structures = reduced_dir / f'{reduced_dir.name}_structures'
reduced_monomers   = reduced_dir / f'{reduced_dir.name}_monomers'

reduced_dir.mkdir(       exist_ok=True)
reduced_monomers.mkdir(  exist_ok=True)
reduced_structures.mkdir(exist_ok=True)

# defining filters
filters = [is_unsolvated]

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    src_mgr = PolymerManager(source_path)

    generate_reduced_structs= src_mgr.logging_wrapper(
        proc_name='Structure reduction', 
        filters=filters
    )(generate_reduced_pdb)

    generate_reduced_structs(
        struct_dir=reduced_structures,
        mono_dir=reduced_monomers,
        DOP=args.DOP,
        max_chain_len=args.max_chain_len,
        flip_term_labels=args.flip_term_labels,
        chain_len_limit=args.chain_len_limit
    )