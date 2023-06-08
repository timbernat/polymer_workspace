'''Functions employed in block-modular scripts for Polymer RCT and simulation workflow'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)
main_logger = logging.getLogger(__name__)

from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler
loggers = [main_logger, *LOGGERS_MASTER]

# Generic imports
from pathlib import Path

# Resource files
import importlib_resources as impres
import resources
avail_sim_templates = resources.AVAIL_RESOURCES['sim_templates']

# Polymer Imports
from polysaccharide.simulation.records import SimulationParameters
from polysaccharide.simulation.execution import run_simulation

from polysaccharide.polymer.representation import Polymer
from polysaccharide.polymer.management import PolymerManager
from polysaccharide.polymer.filters import is_solvated, is_unsolvated, is_charged, filter_factory_by_attr

# Typing and subclassing
from typing import Any, Callable, Iterable, Optional, TypeAlias
PolymerFunct : TypeAlias = Callable[[Polymer, Any], None] # typing schema for a serially-executable polymer function 

# Cheminformatics
from rdkit import Chem
from rdkit.Chem.rdchem import Mol as RDMol

# Molecular Dynamics
from openff.toolkit import ForceField
from openff.interchange import Interchange
from openff.toolkit.topology import Topology, Molecule
from openff.toolkit.typing.engines.smirnoff.parameters import LibraryChargeHandler

from openmm.unit import angstrom, nanometer

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)
SIM_PARAM_PATH = impres.files(resources.sim_templates)


# Function definitions
def assign_polymer_charges(polymer : Polymer) -> None:
    pass

def simulate_polymer(polymer : Polymer, sim_params : SimulationParameters) -> None:
    '''Run OpenMM simulation according to a set of predefined simulation parameters'''
    interchange = polymer.interchange(
        forcefield_path=sim_params.forcefield_path,
        charge_method=sim_params.charge_method,
        periodic=sim_params.periodic
    )

    sim_folder = polymer.make_sim_dir()
    run_simulation(interchange, sim_params=sim_params, output_folder=sim_folder, output_name=polymer.mol_name)