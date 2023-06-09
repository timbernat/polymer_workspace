'''Computes RDF and property time series data and saving to csvs for plotting and analysis'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path
from openmm.unit import nanosecond

# Polymer Imports
from polysaccharide.polymer.management import PolymerManager
from polysaccharide.polymer.filtering import has_sims, is_solvated, is_unsolvated, filter_factory_by_attr

# Utility function imports
from workflow_functs import perform_prop_analysis

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__
)
parser.add_argument('-src', '--source_name'        , help='The name of the target collection of Polymers', required=True)
parser.add_argument('-n', '--mol_names'            , help='If set, will charge ONLY the molecules with the names specified', action='store', nargs='+')
parser.add_argument('-s', '--solv_type'            , help='Set which solvation type to filter for (options are "solv", "unsolv", or "all", defaults to "all")', choices=('solv', 'unsolv', 'all'), nargs='?', default='all')
parser.add_argument('-t', '--sim_time'             , help='If set, will only analyze trajectories run for this number of nanoseconds', action='store', type=float)
parser.add_argument('-tsi', '--traj_sample_interval', help='How often to sample trajectory frames when loading (equilvalent to "stride" in mdtraj); useful for huge trajectories', action='store', type=int, default=1)

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

src_coll_path = COLL_PATH / args.source_name

## defining mol filters
filters = [has_sims]
if args.mol_names is not None:
    desired_mol = filter_factory_by_attr('base_mol_name', lambda name : name in args.mol_names)
    filters.append(desired_mol)

if args.solv_type == 'unsolv':
    filters.append(is_unsolvated)
elif args.solv_type == 'solv':
    filters.append(is_solvated)
else:
    pass # self-documenting placeholder (doesn;t actually do anything)

has_binary_traj = lambda sim_paths, sim_params : (sim_params.report_to_pdb == False)
sim_dir_filters = [has_binary_traj]

if args.sim_time is not None:
    is_long_sim = lambda sim_paths, sim_params : (sim_params.total_time == args.sim_time*nanosecond)
    sim_dir_filters.append(is_long_sim)

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(src_coll_path)

    perform_prop_analysis = mgr.logging_wrapper(
        proc_name='Trajectory Analysis',
        filters=filters
    )(perform_prop_analysis)

    perform_prop_analysis(
        sim_dir_filters=sim_dir_filters, 
        traj_sample_interval=args.traj_sample_interval
    )
