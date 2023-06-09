'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path

# Resource files
import importlib_resources as impres
import resources
avail_chg_templates = resources.AVAIL_RESOURCES['chg_templates']

# Polymer Imports
from polysaccharide.simulation.records import SimulationParameters

from polysaccharide.polymer.management import PolymerManager
from polysaccharide.polymer.filtering import is_unsolvated, is_charged, is_base

# Utility function imports
from workflow_functs import vacuum_anneal

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
SIM_PARAM_PATH = impres.files(resources.sim_templates)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__ # use script docstring as help description 
)
parser.add_argument('-src', '--source_name' , help='The name of the target collection of Polymers', required=True)
parser.add_argument('-sim', '--sim_params'  , help='Name of the simulation parameters preset file to load for simulation', default='vacuum_anneal.json')
parser.add_argument('-r', '--num_replicates', help='Number of total conformers to generate (not counting the original)', type=int, default=4)

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

## defining paths
source_path = COLL_PATH / args.source_name

sim_param_path = SIM_PARAM_PATH / args.sim_params
if not sim_param_path.suffix:
    sim_param_path = sim_param_path.with_name(f'{sim_param_path.stem}.json') # ensure charge params path has correct extension

## defining mol filters
filters = [is_unsolvated, is_charged, is_base]

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(source_path)
    sim_params = SimulationParameters.from_file(sim_param_path)

    vacuum_anneal = mgr.logging_wrapper(
        proc_name=f'Vacuum-anneal conformer generation',
        filters=filters
    )(vacuum_anneal)

    vacuum_anneal(sim_params, num_new_confs=args.num_replicates)