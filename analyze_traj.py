'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
import pandas as pd
from pathlib import Path
from openmm.unit import nanometer

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

main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER] # loggers from all modules which produce logging output

# ------------------------------------------------------------------------------

@mgr.logging_wrapper(loggers, proc_name='Trajectory Analysis', filters=(has_sims))
def perform_prop_analysis(polymer : Polymer, main_logger : logging.Logger, traj_sample_interval : int=1, plot : bool=False):
    '''Analyze trajectories to obtain polymer property data'''
    # aqcuire files for all information
    for sim_folder, sim_paths_file in polymer.simulation_paths.items():
        sim_paths = SimulationPaths.from_file(sim_paths_file)
        sim_params = SimulationParameters.from_file(sim_paths.sim_params)
        state_data = pd.read_csv(sim_paths.state_data)

        if sim_paths.trajectory.suffix == '.dcd': # only attempt to load compressed binary trajectories
            main_logger.info(f'Found DCD trajectory {sim_paths.trajectory}')
            traj = trajectory.load_traj(sim_paths.trajectory, topo_path=polymer.structure_file, sample_interval=traj_sample_interval, remove_solvent=True)

            # save and plot RDF data
            main_logger.info('Calculating pairwise radial distribution functions')
            rdf_dataframe = trajectory.acquire_rdfs(traj, max_rad=1.0*nanometer)
            rdf_save_path = sim_folder/'rdfs.csv'
            sim_paths.spatial_data = rdf_save_path
            rdf_dataframe.to_csv(rdf_save_path, index=False)

            if plot:
                main_logger.info('Plotting pairwise radial distribution functions')
                rdf_fig, rdf_ax = plotprops.plot_rdfs(rdf_dataframe, scale=15.0)
                rdf_fig.suptitle(f'Pairwise Radial Distribution Functions - {polymer.mol_name}')
                rdf_fig.savefig(sim_folder/f'RDFs.png', bbox_inches='tight')
            plt.close()

            # save and plot property data
            main_logger.info('Calculating polymer shape properties')
            prop_dataframe = trajectory.acquire_time_props(traj, time_points=sim_params.time_points[::traj_sample_interval]) 
            prop_save_path = sim_folder/'time_series.csv'
            sim_paths.time_data = prop_save_path
            prop_dataframe.to_csv(prop_save_path, index=False)

            if plot:
                main_logger.info('Plotting polymer shape properties')
                prop_fig, prop_ax = plotprops.plot_time_props(prop_dataframe, scale=18.0)
                prop_fig.suptitle(f'Polymer Shape Properties - {polymer.mol_name}')
                prop_fig.savefig(sim_folder/f'shape_props.png', bbox_inches='tight')
            plt.close()
            
            sim_paths.to_file(sim_paths_file) # update references to analyzed data files in path file
            main_logger.info(f'Successfully exported trajectory analysis data{" and plots" if plot else ""}')

perform_prop_analysis(main_logger)
