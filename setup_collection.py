# Custom Imports
from polymer_utils.general import timestamp_now
from polymer_utils.logutils import MultiStreamFileHandler
from polymer_utils.filetree import startfile

from polymer_utils.representation import PolymerManager
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
    reset      = False
    purge_sims = False
    purge_logs = False

    poly_source_path = COMPAT_PDB_PATH / 'water_soluble_Colina'
    solv_template    = RESOURCE_PATH/'inp_templates'/'solv_polymer_template_box.inp'
    desired_solvents = (WATER_TIP3P,) # (None,)
    exclusion = 1.0*nanometer

    # Define derived paths and create manager
    collection_path  = COLL_PATH / poly_source_path.name
    structure_path   = poly_source_path / f'{poly_source_path.name}_structures'
    monomer_path     = poly_source_path / f'{poly_source_path.name}_monomers'

    mgr = PolymerManager(collection_path)

    # Perform manager setup / purge actions
    if purge_logs: # NOTE : must be done BEFORE log FileHandler is created, as this will destroy it's output as well
        mgr.purge_logs(really=True)

    creation_logger = logging.getLogger('polymer_setup')
    logfile_path = mgr.log_dir/f'Setup_{timestamp_now()}.log'

    with MultiStreamFileHandler(logfile_path, loggers=[creation_logger, polylogger], proc_name=f'Creation of collection "{mgr.collection_dir.name}"'):
        if reset:
            mgr.purge_collection(really=True, purge_logs=False) # Explicitly DON'T purge logs here (will be done prior to entering log loop)

        if purge_sims:
            mgr.purge_sims(really=True)

        if not mgr.polymers: # will be empty if not yet instantiated or if reset prior
            mgr.populate_collection(struct_dir=structure_path, monomer_dir=monomer_path)
            mgr.solvate_collection(desired_solvents, template_path=solv_template, exclusion=exclusion)
    startfile(mgr.collection_dir)