'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
from pathlib import Path
from openmm.unit import nanometer

# Logging
import logging
logging.basicConfig(level=logging.INFO)
from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler

# Polymer Imports
from polysaccharide.representation import PolymerManager
from polysaccharide.solvation.solvents import WATER_TIP3P

# Static Paths
RESOURCE_PATH = Path('resources')
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')

# ------------------------------------------------------------------------------

# Set parameters here
small_coll = COLL_PATH / 'water_soluble_small'
new_coll = COLL_PATH / 'water_soluble_large'
structure_dir = COMPAT_PDB_PATH / 'water_soluble_polymers' / 'water_soluble_polymers_structures' 

# Solvent parameters
solv_template    = RESOURCE_PATH/'inp_templates'/'solv_polymer_template_box.inp'
desired_solvents = (WATER_TIP3P,)
exclusion = 1*nanometer

if __name__ == '__main__':
    new_coll.mkdir(exist_ok=True)
    mgr_large = PolymerManager(new_coll)
    mgr_small = PolymerManager(small_coll)

    logger = logging.getLogger(__name__)
    loggers = [logger, *LOGGERS_MASTER]

        # load small collection and clone into new folder
    with ProcessLogHandler(filedir=mgr_small.log_dir, loggers=loggers, proc_name='Charge transfer'):
        for polymer in mgr_small.polymers_list:
            if polymer.charges and polymer.has_monomer_data_charged:
                polymer.clone(
                    dest_dir=new_coll,
                    clone_name=polymer.base_mol_name, # exclude solvent (will need to resolvate with new structure later)
                    clone_solvent=False,
                    clone_structures=False,
                    clone_monomers=True, # keep only charged monomer information
                    clone_ff=False,
                    clone_charges=False,
                    clone_sims= False
                )

        # copy large structures over to clones, initialize as PolymerManager
    with ProcessLogHandler(filedir=mgr_large.log_dir, loggers=loggers, proc_name='Large structure alignment'):
        mgr_large.update_collection() # update collection to reflect newly created enlarged copies
        for polymer in mgr_large.polymers_list: # copy large structuers over
            polymer.populate_pdb(structure_dir)
        mgr_large.solvate_collection(solvents=desired_solvents, template_path=solv_template, exclusion=exclusion)