'''For transferring residue-averaged charges from reductions to full-sized collection'''

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
from openmm.unit import nanometer

# Resource files
import importlib_resources as impres
import resources
avail_chg_templates = resources.AVAIL_RESOURCES['chg_templates']

# Polymer Imports
from polysaccharide.solvation.solvent import Solvent
from polysaccharide.solvation.solvents import WATER_TIP3P
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import is_charged, has_monomers_chgd, filter_factory_by_attr

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=__doc__ # use script docstring as help description 
)
parser.add_argument('-src', '--source_name', help='The name of the target collection of reduced Polymers', required=True)
parser.add_argument('-out', '--output_name', help='The name of the output collection of Polymers', required=True)
parser.add_argument('-pdb', '--pdb_source' , help='Directory containing the full-sized, structure source files', required=True)
parser.add_argument('-s', '--solvate'      , help='If set, will also solvate all full-sized structures', action='store_true')

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

# defining paths
source_path = COLL_PATH / args.source_name
output_path = COLL_PATH / args.output_name
output_path.mkdir(exist_ok=True)

poly_source_path = COMPAT_PDB_PATH / args.pdb_source
structure_dir  = poly_source_path / f'{poly_source_path.name}_structures'
monomer_dir    = poly_source_path / f'{poly_source_path.name}_monomers'

# defining filters
filters = [is_charged, has_monomers_chgd]

## Hard-coded parameters for now
solv_template    = RESOURCE_PATH/'inp_templates'/'solv_polymer_template_box.inp'
desired_solvents = (WATER_TIP3P,)

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    # load small collection and clone into new folder (over OLD collection)
    mgr_small = PolymerManager(source_path)
    @mgr_small.logging_wrapper(loggers, proc_name='Charge Transfer', filters=filters)
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

    output_path.mkdir(exist_ok=True)
    transfer_charges(output_path)

    # copy large structures over to clones (over NEW collection), initialize as PolymerManager
    mgr_large = PolymerManager(output_path) # NOTE : load after cloning to ensure collection is updated
    @mgr_large.logging_wrapper(loggers, proc_name='Large structure rectification')
    def assign_large_structures(polymer : Polymer, struct_dir : Path, solvate : bool, desired_solvents : tuple[Solvent], template_path : Path) -> None:
        polymer.populate_pdb(struct_dir) # TOSELF : can't use mgr.populate_collection because this creates a new Polmyer (default None value for charged monomer path overwrite reference)
        if solvate:
            polymer.solvate(desired_solvents, template_path=template_path)

    assign_large_structures(structure_dir, solvate=args.solvate, desired_solvents=desired_solvents, template_path=solv_template)