'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
import pandas as pd
from pathlib import Path
from openmm.unit import nanometer, nanosecond

# Logging
import logging
logging.basicConfig(level=logging.INFO)
import matplotlib.pyplot as plt

from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler

main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
from polysaccharide.representation import Polymer, PolymerManager, has_sims
from polysaccharide.solvation.solvents import WATER_TIP3P
from polysaccharide.simulation import SimulationParameters, SimulationPaths
from polysaccharide.analysis import trajectory, plotprops

# Static Paths
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
COLL_PATH = Path('Collections')

RESOURCE_PATH = Path('resources')
CHG_PARAM_PATH = RESOURCE_PATH / 'chg_templates'
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

src_coll_path = COLL_PATH / 'water_soluble_large'
mgr = PolymerManager(src_coll_path)
desired_solvents = (WATER_TIP3P,)

is_long_sim = lambda sim_paths, sim_params : (sim_params.total_time == 100*nanosecond)
has_binary_traj = lambda sim_paths, sim_params : (sim_params.report_to_pdb == False)

# ------------------------------------------------------------------------------

@mgr.logging_wrapper(loggers, proc_name='Trajectory Analysis', filters=(has_sims,))
def perform_prop_analysis(polymer : Polymer, main_logger : logging.Logger, traj_sample_interval : int=1) -> None:
    '''Analyze trajectories to obtain polymer property data'''
    # aqcuire files for all information
    for sim_dir, (sim_paths, sim_params) in polymer.filter_sim_dirs(conditions=(is_long_sim, has_binary_traj)).items():
        main_logger.info(f'Found trajectory {sim_paths.trajectory}')
        state_data = pd.read_csv(sim_paths.state_data)
        traj = trajectory.load_traj(sim_paths.trajectory, topo_path=polymer.structure_file, sample_interval=traj_sample_interval, remove_solvent=True)

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
