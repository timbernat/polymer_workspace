'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Logging
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

from polysaccharide import LOGGERS_MASTER
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
from polysaccharide.solvation.solvents import WATER_TIP3P
from polysaccharide.representation import Polymer, PolymerManager, filter_factory_by_attr
from polysaccharide.simulation.records import SimulationParameters

# Static Paths
COLL_PATH = Path('Collections')
RESOURCE_PATH = Path('resources')
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

src_coll_path = COLL_PATH / 'simple_polymers'
# src_coll_path = COLL_PATH / 'water_soluble_large'
mgr = PolymerManager(src_coll_path)

desired_solvents = (WATER_TIP3P,)
desired_mols = ('polyvinylchloride_solv_water',)
sim_param_paths = [
    SIM_PARAM_PATH / 'pilot_sim_ABE_avg_dcd copy.json'
]
# sim_param_paths = [
#     SIM_PARAM_PATH / 'long_sim_ABE_avg_dcd.json',
#     SIM_PARAM_PATH / 'long_sim_espaloma_dcd.json',
# ]
# sim_param_paths = [
#     SIM_PARAM_PATH / 'standard_sim_ABE_avg_dcd.json',
#     SIM_PARAM_PATH / 'standard_sim_espaloma_dcd.json',
# ]

# ------------------------------------------------------------------------------

# BEGIN CHARGING / SIM LOOP - Perform charge averaging on all target molecules which don't already have averaged LCs; Load forcefield for those which already do 
solvated = filter_factory_by_attr(attr_name='solvent', condition=lambda solv : solv in desired_solvents)
charged  = filter_factory_by_attr(attr_name='charges')
selected = filter_factory_by_attr(attr_name='mol_name', condition=lambda name : name in desired_mols)

for sim_param_path in sim_param_paths:
    sim_params = SimulationParameters.from_file(sim_param_path)
    proc_name = f'Simulation {sim_params.charge_method}'

    @mgr.logging_wrapper(loggers, proc_name=proc_name, filters=(solvated, charged, selected))
    def simulate(polymer : Polymer, sim_params : SimulationParameters) -> None:
        '''Run single NPT-ensemble simulation'''
        polymer.run_simulation(sim_params, ensemble='NPT')

    simulate(sim_params)