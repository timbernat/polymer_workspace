'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
import pandas as pd
from pathlib import Path
from openmm.unit import nanometer

# Logging
import logging
logging.basicConfig(level=logging.INFO)

from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler

main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
from polysaccharide.logutils import ProcessLogHandler
from polysaccharide.representation import Polymer, PolymerManager
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

def perform_prop_analysis(polymer : Polymer, main_logger : logging.Logger, traj_sample_interval : int, plot : bool=False):
    '''Analyze trajectories to obtain polymer property data'''
    # aqcuire files for all information
    sim_folder = polymer.newest_sim_dir
    main_logger.info(f'Acquiring trajectory data from {sim_folder}')
    sim_paths = SimulationPaths.from_file(next(sim_folder.glob('*_paths.json')))

    sim_params = SimulationParameters.from_file(sim_paths.sim_params)
    state_data = pd.read_csv(sim_paths.state_data)
    traj = trajectory.load_traj(sim_paths.trajectory, topo_path=polymer.structure_file, sample_interval=traj_sample_interval, remove_solvent=True)

    # save and plot RDF data
    main_logger.info('Determining pairwise radial distribution functions')
    rdf_dataframe = trajectory.acquire_rdfs(traj, max_rad=1.0*nanometer)
    rdf_dataframe.to_csv(sim_folder/'rdfs.csv')
    if plot:
        rdf_fig, rdf_ax = plotprops.plot_rdfs(rdf_dataframe, scale=15.0)
        rdf_fig.suptitle(f'Pairwise Radial Distribution Functions - {polymer.mol_name}')
        rdf_fig.savefig(sim_folder/f'RDFs.png', bbox_inches='tight')
    plt.close()

    # save and plot property data
    main_logger.info('Determining polymer shape properties')
    prop_dataframe = trajectory.acquire_time_props(traj, properties=analysis.polyprops.DEFAULT_PROPS, time_points=sim_params.time_points[::traj_sample_interval]) 
    prop_dataframe.to_csv(sim_folder/'time_series.csv')
    if plot:
        prop_fig, prop_ax = plotprops.plot_time_props(prop_dataframe, scale=18.0)
        prop_fig.suptitle(f'Polymer Shape Properties - {polymer.mol_name}')
        prop_fig.savefig(sim_folder/f'shape_props.png', bbox_inches='tight')
    plt.close()
    
    main_logger.info('Successfully exported trajectory analysis data and plots')

# ------------------------------------------------------------------------------

src_coll_path = COLL_PATH / 'water_soluble_large'
mgr = PolymerManager(src_coll_path)

desired_solvents = (WATER_TIP3P,)
sim_param_paths = [
    SIM_PARAM_PATH / 'standard_sim_ABE_avg.json',
    SIM_PARAM_PATH / 'standard_sim_espaloma.json',
]

# ------------------------------------------------------------------------------

# BEGIN CHARGING / SIM LOOP - Perform charge averaging on all target molecules which don't already have averaged LCs; Load forcefield for those which already do 
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER] # loggers from all modules which produce logging output

sample_dirs = [
    polymer 
        for polymer in mgr.polymers_list
            if (polymer.solvent in desired_solvents) and polymer.charges
]

for sim_param_path in sim_param_paths:
    sim_params = SimulationParameters.from_file(sim_param_path)
    proc_name = f'Simulation {sim_params.charge_method} and analysis'

    with ProcessLogHandler(filedir=mgr.log_dir, loggers=loggers, proc_name=proc_name, timestamp=True) as msf_handler:
        for i, polymer in enumerate(sample_dirs):
            main_logger.info(f'Current molecule: "{polymer.mol_name}" ({i + 1}/{len(sample_dirs)})') # +1 converts to more human-readable 1-index for step count
            with msf_handler.subhandler(filedir=polymer.logs, loggers=loggers, proc_name=proc_name, timestamp=True) as subhandler: # also log actions to individual Polymers
                sim_folder = polymer.run_simulation_NPT(sim_params)
                perform_prop_analysis(polymer, main_logger)