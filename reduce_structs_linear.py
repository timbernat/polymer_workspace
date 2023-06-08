'''Creates directory of reduced-chain-length structure PDBs and monomer JSONs from a collection of linear Polymers'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)
main_logger = logging.getLogger(__name__)

from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler
loggers = [main_logger, *LOGGERS_MASTER]

# Generic imports
import argparse
from pathlib import Path
from shutil import copyfile
from typing import Iterable

# Resource files
import importlib_resources as impres
import resources
avail_chg_templates = resources.AVAIL_RESOURCES['chg_templates']

# Polymer Imports
from polysaccharide.polymer.representation import Polymer
from polysaccharide.polymer.management import PolymerManager
from polysaccharide.polymer.filters import is_unsolvated

from polysaccharide.polymer.monomer import estimate_max_DOP, estimate_chain_len
from polysaccharide.polymer.building import build_linear_polymer
from polysaccharide.polymer.exceptions import ExcessiveChainLengthError

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__ # use script docstring as help description 
)
parser.add_argument('-src', '--source_name'    , help='The name of the target collection of Polymers', required=True)
parser.add_argument('-out', '--output_name'    , help='The name of the directory to output generated structures to (within the set PDB folder)')
parser.add_argument('-N', '--max_chain_len'    , help='Maximum number of atoms in any of the reduced chain generated. If this is specified, CANNOT specify DOP', type=int)
parser.add_argument('-D', '--DOP'              , help='The number of monomer units to include in the generated reductions.  If this is specified, CANNOT specify max_chain_len', type=int)
parser.add_argument('-lim', '--chain_len_limit', help='The maximum allowable size for a chain to be built to; any chains attempted to be built larger than this limit will raise an error', type=int, default=300)
parser.add_argument('-f', '--flip_term_labels' , help='Names of the chains on which to reverse the order of head/tail terminal group labels (only works for linear homopolymers!)', action='store', nargs='+', default=tuple())

args = parser.parse_args()
if args.output_name is None:
    args.output_name = f'{args.source_name}_reduced' # can't just set as generic default, since this value depends on the source input

# Arg processing
# ------------------------------------------------------------------------------
## Guarantee that exactly one of the mutually exclusive chain length specification is provided
if not (args.max_chain_len or args.DOP):
    raise ValueError('Must provide EITHER a maximum chain length OR a degree of polymerization (provided neither)')

if args.max_chain_len and args.DOP:
    raise ValueError('Must provide EITHER a maximum chain length OR a degree of polymerization (provided both)')

# defining paths
source_path = COLL_PATH / args.source_name

reduced_dir = COMPAT_PDB_PATH / args.output_name # CRITICAL : Note that this is a PDB STRUCTURES path, NOT a COLLECTION path!!
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

    @src_mgr.logging_wrapper(loggers, proc_name='Structure reduction', filters=filters)
    def generate_reduced_pdbs(polymer : Polymer, flip_term_labels : Iterable[str], chain_len_limit : int):
        monomer_smarts = polymer.monomer_info.monomers # create copy to avoid popping from original
        if args.DOP: # NOTE : this only works as intended because of the exclusivity check during arg processing
            DOP = args.DOP # TODO : inject as argument rather than calling global
            max_chain_len = estimate_chain_len(monomer_smarts, DOP)
        if args.max_chain_len:
            max_chain_len = args.max_chain_len
            DOP = estimate_max_DOP(monomer_smarts, max_chain_len)
        
        if max_chain_len > chain_len_limit:
            raise ExcessiveChainLengthError(f'Cannot create reduction with over {chain_len_limit} atoms (requested {max_chain_len})')
        
        chain = build_linear_polymer(monomer_smarts, DOP=DOP, reverse_term_labels=(polymer.mol_name in flip_term_labels))
        chain.save(str(reduced_structures/f'{polymer.mol_name}.pdb'), overwrite=True)
        copyfile(polymer.monomer_file_uncharged, reduced_monomers/f'{polymer.mol_name}.json')

    generate_reduced_pdbs(flip_term_labels=args.flip_term_labels, chain_len_limit=args.chain_len_limit)