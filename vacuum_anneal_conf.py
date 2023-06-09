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

from polysaccharide.polymer.management import PolymerManager, MolFilterBuffer
from polysaccharide.polymer.filtering import is_base

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
MolFilterBuffer.argparse_inject(parser)

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

## defining paths
source_path = COLL_PATH / args.source_name

sim_param_path = SIM_PARAM_PATH / args.sim_params
if not sim_param_path.suffix:
    sim_param_path = sim_param_path.with_name(f'{sim_param_path.stem}.json') # ensure charge params path has correct extension

## defining mol filters
molbuf = MolFilterBuffer.from_argparse(args)
molbuf.solvent = False # force preference for unsolvated molecules (VACUUM anneal)
molbuf.charges = True  # force preference for charged molecules (can;t run sim otherwise)
mol_filters = molbuf.filters
mol_filters.append(is_base) # also assert that only base molecules are annealed

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(source_path)
    sim_params = SimulationParameters.from_file(sim_param_path)

    vacuum_anneal = mgr.logging_wrapper(
        proc_name=f'Vacuum-anneal conformer generation',
        filters=mol_filters
    )(vacuum_anneal)

    vacuum_anneal(sim_params, num_new_confs=args.num_replicates)