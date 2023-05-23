'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

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
from polysaccharide.representation import is_solvated, is_unsolvated, is_uncharged, filter_factory_by_attr
from polysaccharide.charging.application import ChargingParameters

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description='Assigns partial charges to the Polymers in a collection, with optional filtering'
)
parser.add_argument('-src', '--source_name' , help='The name of the target collection of Polymers', required=True)
parser.add_argument('-cp', '--charge_params', help=f'Name of the charging parameters preset file to load for charging (available files are {", ".join(avail_chg_templates)})', required=True)
parser.add_argument('-n', '--mol_names'     , help='If set, will charge ONLY the molecules with the names specified', action='store', nargs='+')
parser.add_argument('-s', '--solv_type'     , help='Set which solvation type to filter for (options are "solv", "unsolv", or "all", defaults to "unsolv")', choices=('solv', 'unsolv', 'all'), nargs='?', const='unsolv')

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

## defining paths
src_coll_path = COLL_PATH / args.source_name

chg_params_path = impres.files(resources.chg_templates) / args.charge_params
if not chg_params_path.suffix:
    chg_params_path = chg_params_path.with_name(f'{chg_params_path.stem}.json') # ensure charge params path has correct extension
chg_params = ChargingParameters.from_file(chg_params_path)

## defining mol filters
filters = [is_uncharged]
if args.mol_names is not None:
    desired_mol = filter_factory_by_attr('base_mol_name', lambda name : name in args.mol_names)
    filters.append(desired_mol)

if args.solv_type == 'unsolv':
    filters.append(is_unsolvated)
elif args.solv_type == 'solv':
    filters.append(is_solvated)
else:
    pass # self-documenting placeholder (doesn;t actually do anything)

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(src_coll_path)
    generate_charge_files = mgr.logging_wrapper(
        loggers,
        proc_name=f'Charge assignment',
        filters=filters
    )(Polymer.obtain_partial_charges) 
    
    generate_charge_files(chg_params)
