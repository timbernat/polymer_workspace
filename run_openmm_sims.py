'''Dispatches chosen charged molecules to OpenMM for simulation'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path

# Resource files
import importlib_resources as impres
import resources
avail_sim_templates = ', '.join(
    path.name
        for path in resources.AVAIL_RESOURCES['sim_templates']
)

# Polymer Imports
from polysaccharide.simulation.records import SimulationParameters

from polysaccharide.polymer.representation import Polymer
from polysaccharide.polymer.management import PolymerManager
from polysaccharide.polymer.filtering import is_solvated, is_unsolvated, is_charged, filter_factory_by_attr

# Utility function imports
from workflow_functs import simulate_polymer

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
SIM_PARAM_PATH = impres.files(resources.sim_templates)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__
)
parser.add_argument('-src', '--source_name', help='The name of the target collection of Polymers', required=True)
parser.add_argument('-sp', '--sim_params'  , help=f'Name of the simulation parameters preset file to load for simulation (available files are {avail_sim_templates})', action='store', nargs='+', required=True)
parser.add_argument('-n', '--mol_names'    , help='If set, will charge ONLY the molecules with the names specified', action='store', nargs='+')
parser.add_argument('-s', '--solv_type'    , help='Set which solvation type to filter for (options are "solv", "unsolv", or "all", defaults to "solv")', choices=('solv', 'unsolv', 'all'), nargs='?', default='solv')

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

## defining paths
src_coll_path = COLL_PATH / args.source_name

sim_param_paths = []
for sim_param_name in args.sim_params:
    sim_param_path = SIM_PARAM_PATH / sim_param_name
    if not sim_param_path.suffix:
        sim_param_path = sim_param_path.with_name(f'{sim_param_path.stem}.json') # ensure charge params path has correct extension
    sim_param_paths.append(sim_param_path)

## defining mol filters
filters = [is_charged]
if args.mol_names is not None:
    desired_mol = filter_factory_by_attr('base_mol_name', lambda name : name in args.mol_names)
    filters.append(desired_mol)

if args.solv_type == 'unsolv':
    filters.append(is_unsolvated)
elif args.solv_type == 'solv':
    filters.append(is_solvated)
else:
    pass # self-documenting placeholder (doesn;t actually do anything)

# ------------------------------------------------------------------------------

# BEGIN CHARGING / SIM LOOP - Perform charge averaging on all target molecules which don't already have averaged LCs; Load forcefield for those which already do 
if __name__ == '__main__':
    mgr = PolymerManager(src_coll_path)
    for sim_param_path in sim_param_paths:
        simulate = mgr.logging_wrapper(
            proc_name=f'Simulation {sim_param_path.stem}',
            filters=filters
        )(simulate_polymer)

        sim_params = SimulationParameters.from_file(sim_param_path)
        simulate(sim_params)