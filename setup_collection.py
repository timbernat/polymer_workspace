# Custom Imports
from polymer_utils.general import timestamp_now
from polymer_utils.logutils import config_mlf_handler
from polymer_utils.filetree import startfile

from polymer_utils.representation import PolymerDirManager
from polymer_utils.representation import LOGGER as polylogger
from polymer_utils.solvation.solvents import WATER_TIP3P

# General Imports
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
                            
from openmm.unit import nanometer

# Static Paths
RESOURCE_PATH = Path('resources')
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs')

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    # Configure paths and parameters
    reset      = True
    purge_sims = False 
    purge_logs = True

    poly_source_path = COMPAT_PDB_PATH / 'testing'
    structure_path   = poly_source_path / 'testing_structures'
    monomer_path     = poly_source_path / 'testing_monomers'
    collection_path  = COLL_PATH / poly_source_path.name
    solv_template    = RESOURCE_PATH/'inp_templates'/'solv_polymer_template_box.inp'

    desired_solvents = (WATER_TIP3P,) # (None,)
    exclusion = 1.0*nanometer

    # Create manager and set up directories, logging all
    mgr = PolymerDirManager(collection_path)

    if purge_logs: # imperative that this be done BEFORE the multi-file logger is created (would clear new log otherwise)
        mgr.purge_logs(really=True)

    creation_logger = logging.getLogger('polymer_setup')
    loggers = [creation_logger, polylogger] 
    creation_log_handler = config_mlf_handler(mgr.log_dir/f'Setup_{timestamp_now()}.log', loggers, writemode='a')

    if reset:
        mgr.purge_collection(really=True) 

    if purge_sims:
        mgr.purge_sims(really=True)

    if not mgr.mol_dirs: # will be empty if not yet instantiated or if reset prior
        mgr.populate_collection(struct_dir=structure_path, monomer_dir=monomer_path)
        mgr.solvate_collection(desired_solvents, template_path=solv_template, exclusion=exclusion)
        creation_logger.info(f'Completed creation of collection "{collection_path.name}"')
        
    creation_log_handler.remove_from_loggers(*loggers)
    startfile(mgr.collection_dir)