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
from polysaccharide.logutils import ProcessLogHandler
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.solvation.solvents import WATER_TIP3P
from polysaccharide.simulation import SimulationParameters, SimulationPaths
from polysaccharide.analysis import trajectory, plotprops, polyprops

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