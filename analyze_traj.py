'''Computes RDF and property time series data and saving to csvs for plotting and analysis'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)
main_logger = logging.getLogger(__name__)

from polysaccharide import LOGGERS_MASTER
loggers = [main_logger, *LOGGERS_MASTER]

# Generic imports
import argparse
from pathlib import Path
import pandas as pd
from openmm.unit import nanometer, nanosecond

# Resource files
import importlib_resources as impres
import resources
avail_sim_templates = resources.AVAIL_RESOURCES['sim_templates']

# Polymer Imports
from polysaccharide.analysis import trajectory
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import has_sims, is_solvated, is_unsolvated, filter_factory_by_attr

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__
)
parser.add_argument('-src', '--source_name', help='The name of the target collection of Polymers', required=True)
parser.add_argument('-n', '--mol_names'    , help='If set, will charge ONLY the molecules with the names specified', action='store', nargs='+')
parser.add_argument('-s', '--solv_type'    , help='Set which solvation type to filter for (options are "solv", "unsolv", or "all", defaults to "all")', choices=('solv', 'unsolv', 'all'), nargs='?', const='all')
parser.add_argument('-t', '--sim_time'     , help='If set, will only analyze trajectories run for this number of nanoseconds', action='store', type=float)

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

    @mgr.logging_wrapper(loggers, proc_name='Trajectory Analysis', filters=filters)
    def perform_prop_analysis(polymer : Polymer, main_logger : logging.Logger, traj_sample_interval : int=1) -> None:
        '''Analyze trajectories to obtain polymer property data'''
        # aqcuire files for all information
        for sim_dir, (sim_paths, sim_params) in polymer.filter_sim_dirs(conditions=sim_dir_filters).items():
            main_logger.info(f'Found trajectory {sim_paths.trajectory}')
            traj = polymer.load_traj(sim_dir)

            # save and plot RDF data
            main_logger.info('Calculating pairwise radial distribution functions')
            rdf_dataframe = trajectory.acquire_rdfs(traj, max_rad=1.0*nanometer)
            rdf_save_path = sim_dir/'rdfs.csv'
            sim_paths.spatial_data = rdf_save_path
            rdf_dataframe.to_csv(rdf_save_path, index=False)

            # save and plot property data
            main_logger.info('Calculating polymer shape properties')
            prop_dataframe = trajectory.acquire_time_props(traj, time_points=sim_params.time_points[::traj_sample_interval]) 
            prop_save_path = sim_dir/'time_series.csv'
            sim_paths.time_data = prop_save_path
            prop_dataframe.to_csv(prop_save_path, index=False)

            sim_paths.to_file(polymer.simulation_paths[sim_dir]) # update references to analyzed data files in path file
            main_logger.info(f'Successfully exported trajectory analysis data')

    perform_prop_analysis(main_logger)
