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
from typing import Iterable

# Resource files
import importlib_resources as impres
import resources
avail_chg_templates = resources.AVAIL_RESOURCES['chg_templates']

# Custom Imports
from polysaccharide import LOGGERS_MASTER
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import is_unsolvated, filter_factory_by_attr
from polysaccharide.solvation import solvents as psolvents
from polysaccharide.solvation.solvent import Solvent
from openmm.unit import nanometer

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)
SOLV_TEMP_PATH = impres.files(resources.inp_templates)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__
)
parser.add_argument('-src', '--source_name', help='The name of the target collection of Polymers', required=True)
parser.add_argument('-n', '--mol_names'    , help='If set, will charge ONLY the molecules with the names specified', action='store', nargs='+')
parser.add_argument('-s' , '--solvents'    , help='Names of all solvent molecule to solvate the target systems in', action='store', nargs='+')
parser.add_argument('-t' , '--template'    , help='Name of the packmol input template file to use for solvation', action='store', default='solv_polymer_template_box.inp')
parser.add_argument('-e' , '--exclusion'   , help='Distance (in nm) between the bounding box of the molecule and the simiulation / solvation box', action='store', type=int, default=1)

args = parser.parse_args()

# Fixed parameters
# ------------------------------------------------------------------------------
exclusion = 1.0*nanometer

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
filters = [is_unsolvated]
if args.mol_names is not None:
    desired_mol = filter_factory_by_attr('base_mol_name', lambda name : name in args.mol_names)
    filters.append(desired_mol)

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(source_path)

    @mgr.logging_wrapper(loggers, proc_name='Solvation', filters=filters)
    def solvate_collection(polymer : Polymer, solvents : Iterable[Solvent], template_path : Path, exclusion : float) -> None:
        polymer.solvate(solvents, template_path=template_path, exclusion=exclusion)

    solvate_collection(solvents=desired_solvents, template_path=solv_template, exclusion=exclusion)