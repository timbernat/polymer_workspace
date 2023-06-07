'''(Re)populate a collection of Polymer directories from a source directory of structural and monomer information'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)
main_logger = logging.getLogger(__name__)

from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler
loggers = [main_logger, *LOGGERS_MASTER]

# Generic imports
import argparse
from pathlib import Path

# Resource files
import importlib_resources as impres
import resources
avail_chg_templates = resources.AVAIL_RESOURCES['chg_templates']

# Custom Imports
from polysaccharide import LOGGERS_MASTER
from polysaccharide.representation import PolymerManager
from polysaccharide.solvation.solvents import WATER_TIP3P
from openmm.unit import nanometer

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__
)
parser.add_argument('-strct', '--structures', help='The name of the directory in the set compatible_pdbs folder containing the target pdb and monomer directories', required=True)
parser.add_argument('-out', '--output_name' , help='The name of the output collection of Polymers')
parser.add_argument('-r'  , '--reset'       , help='If set, will delete the target collection if it already exists'        , action='store_true')
parser.add_argument('-ps' , '--purge_sims'  , help='If set, will delete any extant MD simulations in the target collection', action='store_true')
parser.add_argument('-pl' , '--purge_logs'  , help='If set, will delete any extant log files in the target collection'     , action='store_true')

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

## defining paths
poly_source_path = COMPAT_PDB_PATH / args.structures
structure_path   = poly_source_path / f'{poly_source_path.name}_structures'
monomer_path     = poly_source_path / f'{poly_source_path.name}_monomers'

if args.output_name is None:
    args.output_name = args.structures
collection_path  = COLL_PATH / args.output_name

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(collection_path)

    # Purge logs (if desired), then set up new loggesr
    if args.purge_logs: # NOTE : must be done BEFORE log FileHandler is created, as this will destroy it's output as well
        mgr.purge_logs(really=True)

    # Perform manager file setup / purge actions
    with ProcessLogHandler(filedir=mgr.log_dir, loggers=loggers, proc_name=f'Setup of {mgr.collection_dir.name}', timestamp=True):
        if args.reset:
            mgr.purge_collection(really=True, purge_logs=False) # Explicitly DON'T purge logs here (will be done prior to entering log loop)

        if args.purge_sims:
            mgr.purge_sims(really=True)

        if not mgr.polymers: # will be empty if not yet instantiated or if reset prior
            mgr.populate_collection(struct_dir=structure_path, monomer_dir=monomer_path)