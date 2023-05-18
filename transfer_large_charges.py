'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''
'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
from pathlib import Path
from openmm.unit import nanometer
from pathlib import Path

# Logging
import logging
logging.basicConfig(level=logging.INFO)

from polysaccharide import LOGGERS_MASTER
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
from polysaccharide.solvation.solvent import Solvent
from polysaccharide.solvation.solvents import WATER_TIP3P
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import is_charged, has_monomers_chgd

# Static Paths
COLL_PATH = Path('Collections')
RESOURCE_PATH = Path('resources')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

# Set parameters here
small_coll = COLL_PATH / 'water_soluble_reduced'
new_coll = COLL_PATH / 'water_soluble_large_conf'

src_file_dir = COMPAT_PDB_PATH / 'water_soluble_polymers'
structure_dir = src_file_dir / 'water_soluble_polymers_structures' 
monomer_dir  = src_file_dir  / 'water_soluble_polymers_monomers' 

# Solvent parameters
solv_template    = RESOURCE_PATH/'inp_templates'/'solv_polymer_template_box.inp'
desired_solvents = (WATER_TIP3P,)
exclusion = 1*nanometer

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    # load small collection and clone into new folder
    mgr_small = PolymerManager(small_coll)
    @mgr_small.logging_wrapper(loggers, proc_name='Charge Transfer', filters=(is_charged, has_monomers_chgd))
    def transfer_charges(polymer : Polymer, collection_path : Path) -> None:
        dest_dir = collection_path/polymer.base_mol_name
        dest_dir.mkdir(exist_ok=True)

        polymer.clone(
            dest_dir=dest_dir,
            clone_name=polymer.base_mol_name, 
            clone_solvent=False, # exclude solvent (will need to resolvate with new structure later anyway)
            clone_structures=False,
            clone_monomers=True, # keep only charged monomer information
            clone_ff=False,
            clone_charges=False,
            clone_sims= False
        )

    new_coll.mkdir(exist_ok=True)
    transfer_charges(new_coll)

    # copy large structures over to clones, initialize as PolymerManager
    mgr_large = PolymerManager(new_coll) # NOTE : load after cloning to ensure collection is updated
    @mgr_large.logging_wrapper(loggers, proc_name='Large structure rectification')
    def assign_large_structures(polymer : Polymer, struct_dir : Path, desired_solvents : tuple[Solvent], template_path : Path, exclusion : float) -> None:
        polymer.populate_pdb(struct_dir) # TOSELF : can't use mgr.populate_collection because this creates a new Polmyer (default None value for charged monomer path overwrite reference)
        polymer.solvate_multi(desired_solvents, template_path=template_path, exclusion=exclusion)

    assign_large_structures(structure_dir, desired_solvents, solv_template, exclusion)