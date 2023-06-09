'''(Re)populate a collection of Polymer directories from a source directory of structural and monomer information'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path
from openmm.unit import nanometer

# Resource files
import importlib_resources as impres
import resources
avail_chg_templates = resources.AVAIL_RESOURCES['chg_templates']

# Custom Imports
from polysaccharide.polymer.management import PolymerManager, MolFilterBuffer
from polysaccharide.solvation import solvents as psolvents
from polysaccharide.solvation.solvent import Solvent

# Utility function imports
from workflow_functs import solvate

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
SOLV_TEMP_PATH = impres.files(resources.inp_templates)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__
)
parser.add_argument('-src', '--source_name', help='The name of the target collection of Polymers', required=True)
parser.add_argument('-s', '--solvents'     , help='Names of all solvent molecule to solvate the target systems in' , action='store', nargs='+', default=['WATER_TIP3P'])
parser.add_argument('-t', '--template'     , help='Name of the packmol input template file to use for solvation', action='store', default='solv_polymer_template_box.inp')
parser.add_argument('-e', '--exclusion'    , help='Distance (in nm) between the bounding box of the molecule and the simiulation / solvation box', action='store', type=float, default=1.0)
MolFilterBuffer.argparse_inject(parser)

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

## defining paths
source_path   = COLL_PATH / args.source_name
solv_template = SOLV_TEMP_PATH / args.template

## defining solvents and box volumes
if not args.solvents:
    raise ValueError('Must specify at least 1 solvent')
desired_solvents = [
    getattr(psolvents, solvent_name)
        for solvent_name in args.solvents
]

exclusion = args.exclusion * nanometer # assign units

## defining filters
molbuf = MolFilterBuffer.from_argparse(args)
molbuf.solvent = False # force preference for only unsolvated molecules (don't want to attempt solvation twice)

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(source_path)

    print(molbuf, molbuf.valid_names(mgr))
    solvate_collection = mgr.logging_wrapper(
        proc_name='Solvation',
        filters=molbuf.filters
    )(solvate)

    solvate_collection(solvents=desired_solvents, template_path=solv_template, exclusion=exclusion)