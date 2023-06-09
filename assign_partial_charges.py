'''Assigns partial charges to the Polymers in a collection, with optional filtering by name and solvent type'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path

# Resource files
import importlib_resources as impres
import resources
avail_chg_templates = ', '.join(
    path.name
        for path in resources.AVAIL_RESOURCES['chg_templates']
)

# Polymer Imports
from polysaccharide.charging.application import ChargingParameters
from polysaccharide.polymer.management import PolymerManager, MolFilterBuffer

# Utility function imports
from workflow_functs import assign_polymer_charges

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__
)
parser.add_argument('-src', '--source_name' , help='The name of the target collection of Polymers', required=True)
parser.add_argument('-cp', '--charge_params', help=f'Name of the charging parameters preset file to load for charging (available files are {avail_chg_templates})', required=True)
MolFilterBuffer.argparse_inject(parser)

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
molbuf = MolFilterBuffer.from_argparse(args)
# TOSELF : overwrite / charge status force here?

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(src_coll_path)
    assign_charges = mgr.logging_wrapper(
        proc_name=f'Charge assignment',
        filters=molbuf.filters
    )(assign_polymer_charges) 
    
    assign_charges(chg_params)
