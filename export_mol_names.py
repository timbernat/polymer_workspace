'''Generic template for defining CLI-executable polymer action scripts'''

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

# Resource files
import importlib_resources as impres
import resources
avail_chg_templates = resources.AVAIL_RESOURCES['chg_templates']

# Polymer Imports
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import is_solvated, is_unsolvated, is_charged, is_uncharged, filter_factory_by_attr
from polysaccharide.charging.application import ChargingParameters

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__ # use script docstring as help description 
)
parser.add_argument('-src', '--source_name' , help='The name of the target collection of Polymers', required=True)
parser.add_argument('-out', '--output_file' , help='Name of the text file to output polymer names to', required=True)
parser.add_argument('-n', '--mol_names'     , help='If set, will charge ONLY the molecules with the names specified', action='store', nargs='+')
parser.add_argument('-s', '--solv_type'     , help='Set which solvation type to filter for (options are "solv", "unsolv", or "all", defaults to "all")', choices=('solv', 'unsolv', 'all'), nargs='?', default='all')
parser.add_argument('-c', '--charge_type'   , help='Set which charging status to filter for (options are "chg", "unchg", or "all", defaults to "all")' , choices=('chg', 'unchg', 'all')  , nargs='?', default='all')

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

## defining paths
source_path = COLL_PATH / args.source_name

name_output_path = Path(args.output_file)
if not name_output_path.suffix:
    name_output_path = name_output_path.with_name(f'{name_output_path.stem}.txt') # ensure charge params path has correct extension
name_output_path.touch()

## defining mol filters
filters = []
if args.mol_names is not None:
    desired_mol = filter_factory_by_attr('base_mol_name', lambda name : name in args.mol_names)
    filters.append(desired_mol)

if args.solv_type == 'unsolv':
    filters.append(is_unsolvated)
elif args.solv_type == 'solv':
    filters.append(is_solvated)
else:
    pass # self-documenting placeholder (doesn;t actually do anything)

if args.solv_type == 'unchg':
    filters.append(is_uncharged)
elif args.solv_type == 'chg':
    filters.append(is_charged)
else:
    pass # self-documenting placeholder (doesn;t actually do anything)

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(source_path)
    with name_output_path.open('w') as name_output_file:
        for mol_name in mgr.filtered_by(filters):
            name_output_file.write(f'{mol_name}\n')