'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Logging
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

from polysaccharide import LOGGERS_MASTER
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import is_unsolvated, is_charged, is_base
from polysaccharide.simulation.records import SimulationParameters

# Static Paths
COLL_PATH = Path('Collections')
RESOURCE_PATH = Path('resources')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

src_coll_path = COLL_PATH / 'water_soluble_large_conf'
sim_param_path = SIM_PARAM_PATH / 'vacuum_anneal_sim_dcd.json'
num_confs = 2

# ------------------------------------------------------------------------------

# BEGIN CHARGING / SIM LOOP - Perform charge averaging on all target molecules which don't already have averaged LCs; Load forcefield for those which already do 
if __name__ == '__main__':
    mgr = PolymerManager(src_coll_path)
    sim_params = SimulationParameters.from_file(sim_param_path)

    @mgr.logging_wrapper(loggers, proc_name=f'Vacuum-anneal conformer generation', filters=(is_unsolvated, is_charged, is_base))
    def vacuum_anneal(polymer : Polymer, main_logger : logging.Logger, sim_params : SimulationParameters, num_confs : int, snapshot_idx : int=-1) -> None:
        '''Run quick vacuum NVT sim at high T'''
        for i in range(num_confs):
            conf_clone = polymer.clone(
                clone_name=f'{polymer.base_mol_name}_conf_{i + 1}',
                clone_solvent=True,
                clone_structures=True,
                clone_monomers=True,
                clone_ff=True,
                clone_charges=True,
                clone_sims=False
            )

            polymer.run_simulation(sim_params, ensemble='NVT', periodic=False)
            
            main_logger.info('Extracting final conformation from simulation')
            traj = polymer.load_traj(polymer.newest_sim_dir)
            new_conf = traj[snapshot_idx]
            main_logger.info('Applying new conformation to clone')
            new_conf.save(conf_clone.structure_file) # overwrite the clone's structure with the new conformer

    vacuum_anneal(main_logger, sim_params, num_confs=num_confs)