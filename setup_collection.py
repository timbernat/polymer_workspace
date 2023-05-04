'''For initializing a managed collection of Polymers from directories of structure and monomer files'''

# Custom Imports
from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler
from polysaccharide.filetree import startfile
from polysaccharide.representation import PolymerManager
from polysaccharide.solvation.solvents import WATER_TIP3P

# Generic Imports
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
                            
from openmm.unit import nanometer

# Static Paths
RESOURCE_PATH = Path('resources')
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')

# ------------------------------------------------------------------------------

# Set parameters here
reset      = True
purge_sims = True 
purge_logs = True

# poly_source_path = COMPAT_PDB_PATH / 'simple_polymers'
# poly_source_path = COMPAT_PDB_PATH / 'water_soluble_polymers'
poly_source_path = COMPAT_PDB_PATH / 'water_soluble_small'
solv_template    = RESOURCE_PATH/'inp_templates'/'solv_polymer_template_box.inp'
desired_solvents = (WATER_TIP3P,) # (None,)
exclusion = 1.0*nanometer

if __name__ == '__main__':
    # Define derived paths and create manager
    collection_path  = COLL_PATH / poly_source_path.name
    structure_path   = poly_source_path / f'{poly_source_path.name}_structures'
    monomer_path     = poly_source_path / f'{poly_source_path.name}_monomers'

    mgr = PolymerManager(collection_path)

    # Perform manager setup / purge actions
    if purge_logs: # NOTE : must be done BEFORE log FileHandler is created, as this will destroy it's output as well
        mgr.purge_logs(really=True)

    creation_logger = logging.getLogger('polymer_setup')
    loggers = [creation_logger, *LOGGERS_MASTER]

    with ProcessLogHandler(filedir=mgr.log_dir, loggers=loggers, proc_name=f'Setup of {mgr.collection_dir.name}', timestamp=True):
        if reset:
            mgr.purge_collection(really=True, purge_logs=False) # Explicitly DON'T purge logs here (will be done prior to entering log loop)

        if purge_sims:
            mgr.purge_sims(really=True)

        if not mgr.polymers: # will be empty if not yet instantiated or if reset prior
            mgr.populate_collection(struct_dir=structure_path, monomer_dir=monomer_path)
            mgr.solvate_collection(desired_solvents, template_path=solv_template, exclusion=exclusion)