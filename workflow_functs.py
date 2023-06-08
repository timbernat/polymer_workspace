'''Functions employed in block-modular scripts for Polymer RCT and simulation workflow'''

# Logging
import logging
from logging import Logger
logging.basicConfig(level=logging.INFO, force=True)

from polysaccharide import LOGGERS_MASTER

# Generic imports
from pathlib import Path
from shutil import copyfile

# Polymer Imports
from polysaccharide.solvation.solvent import Solvent

from polysaccharide.charging.application import ChargingParameters, CHARGER_REGISTRY
from polysaccharide.charging.averaging import get_averaged_residue_charges, AveragingCharger

from polysaccharide.simulation.records import SimulationParameters
from polysaccharide.simulation.execution import run_simulation

from polysaccharide.polymer.representation import Polymer, SimDirFilter
from polysaccharide.polymer.management import PolymerManager, PolymerFunction
from polysaccharide.polymer.filters import is_solvated, is_unsolvated, is_charged, filter_factory_by_attr
from polysaccharide.polymer.monomer import estimate_max_DOP, estimate_chain_len
from polysaccharide.polymer.building import build_linear_polymer
from polysaccharide.polymer.exceptions import ExcessiveChainLengthError

from polysaccharide.analysis import trajectory

# Typing and subclassing
from typing import Any, Callable, Iterable, Optional, TypeAlias, Union

# Cheminformatics
from rdkit import Chem
from rdkit.Chem.rdchem import Mol as RDMol

# Molecular Dynamics
from openff.toolkit import ForceField
from openff.interchange import Interchange
from openff.toolkit.topology import Topology, Molecule
from openff.toolkit.typing.engines.smirnoff.parameters import LibraryChargeHandler

from openmm.unit import angstrom, nanometer

# Resource files
import importlib_resources as impres
import resources

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
SIM_PARAM_PATH = impres.files(resources.sim_templates)


# Function definitions
def PLACEHOLDER(polymer : Polymer, poly_logger : Logger) -> None:
    '''Placeholder template function - for self-documentation AND for script_template.py'''
    pass

def solvate(polymer : Polymer, poly_logger : Logger, solvents : Union[Solvent, Iterable[Solvent]], template_path : Path, exclusion : float=None) -> None:
    '''Fill a box around a polymer with solvent'''
    polymer.solvate(solvents, template_path=template_path, exclusion=exclusion)

def generate_reduced_pdb(polymer : Polymer, poly_logger : Logger, struct_dir : Path, mono_dir : Path, DOP : Optional[int], max_chain_len : Optional[int], flip_term_labels : Iterable[str], chain_len_limit : int=300):
    '''Builds new PDB structures of the desired size from the monomers of an existing Polymer'''
    monomer_smarts = polymer.monomer_info.monomers # create copy to avoid popping from original
    if not (max_chain_len or DOP):
        raise ValueError('Must provide EITHER a maximum chain length OR a degree of polymerization (provided neither)')

    if max_chain_len and DOP:
        raise ValueError('Must provide EITHER a maximum chain length OR a degree of polymerization (provided both)')

    if DOP: # NOTE : this only works as intended because of the exclusivity check during arg processing
        max_chain_len = estimate_chain_len(monomer_smarts, DOP)

    if max_chain_len:
        DOP = estimate_max_DOP(monomer_smarts, max_chain_len)
    
    if max_chain_len > chain_len_limit:
        raise ExcessiveChainLengthError(f'Cannot create reduction with over {chain_len_limit} atoms (requested {max_chain_len})')
    
    chain = build_linear_polymer(monomer_smarts, DOP=DOP, reverse_term_labels=(polymer.mol_name in flip_term_labels))
    chain.save(str(struct_dir/f'{polymer.mol_name}.pdb'), overwrite=True)
    copyfile(polymer.monomer_file_uncharged, mono_dir/f'{polymer.mol_name}.json')

def assign_polymer_charges(polymer : Polymer, poly_logger : Logger, chg_params : ChargingParameters) -> None:
    '''Ensure a Polymer has all partial charge sets'''
    assert(polymer.has_monomer_info)

    # 1) ENSURING CHARGES AND RELATED FILES FOR ALL BASE CHARGING METHODS EXIST
    for chg_method in chg_params.charge_methods:
        chgr = CHARGER_REGISTRY[chg_method]()
        polymer.assert_charges_for(chgr, return_cmol=False)
        poly_logger.info('') # add gaps between charge method for breathing room
        
    # 2) GENERATE RESIDUE-AVERAGED LIBRARY CHARGES FROM THE CHARGE SET OF CHOICE
    if polymer.has_monomer_info_charged: # load precomputed residue charges if existing set is found
        poly_logger.info(f'Found charged JSON, loading averaged charges for {polymer.mol_name} residues from file')
        residue_charges = polymer.monomer_info_charged.charges
    else: # otherwise, recompute them from the charge set of choice
        poly_logger.info(f'Averaging {chg_params.averaging_charge_method} charges over {polymer.mol_name} residues')
        residue_charges = get_averaged_residue_charges(
            cmol=polymer.charged_offmol(chg_params.averaging_charge_method), # TODO : add check to ensure this ISN'T and averaging method
            monomer_info=polymer.monomer_info_uncharged
        )
    
    if (not polymer.has_monomer_info_charged) or chg_params.overwrite_chg_mono: # cast charges residues to file if none exist or is explicitly overwriting
        polymer.create_charged_monomer_file(residue_charges) # TOSELF : this is separate from the above clauses as it might be called regardless of existing charges during overwrite

    avg_chgr = AveragingCharger() # generate precomputed charge set for full molecule from residue library charges
    avg_chgr.set_residue_charges(residue_charges)
    polymer.assert_charges_for(avg_chgr, return_cmol=False)

def retrieve_monomers(polymer : Polymer, poly_logger : Logger, charged_mgr : PolymerManager) -> None:
    '''Copies charged monomer files from a corresponding Polymer in another collection'''
    counterpart = charged_mgr.polymers[polymer.mol_name]
    assert(polymer.solvent == counterpart.solvent)

    if not polymer.has_monomer_info_uncharged:
        counterpart.transfer_file_attr('monomer_file_uncharged', polymer)
    
    if not polymer.has_monomer_info_charged:
        counterpart.transfer_file_attr('monomer_file_charged', polymer)

def simulate_polymer(polymer : Polymer, poly_logger : Logger, sim_params : SimulationParameters) -> None:
    '''Run OpenMM simulation according to a set of predefined simulation parameters'''
    interchange = polymer.interchange(
        forcefield_path=sim_params.forcefield_path,
        charge_method=sim_params.charge_method,
        periodic=sim_params.periodic
    )

    sim_folder = polymer.make_sim_dir()
    run_simulation(interchange, sim_params=sim_params, output_folder=sim_folder, output_name=polymer.mol_name)

def vacuum_anneal(polymer : Polymer, poly_logger : logging.Logger, sim_params : SimulationParameters, num_new_confs : int, snapshot_idx : int=-1) -> None:
    '''Run quick vacuum NVT sim at high T'''
    for i in range(num_new_confs):
        conf_clone = polymer.clone(
            clone_name=f'{polymer.base_mol_name}_conf_{i + 1}',
            clone_solvent=True,
            clone_structures=True,
            clone_monomers=True,
            clone_charges=True,
            clone_sims=False
        )
        simulate_polymer(polymer, sim_params)
        
        poly_logger.info('Extracting final conformation from simulation')
        traj = polymer.load_traj(polymer.newest_sim_dir)
        new_conf = traj[snapshot_idx]
        poly_logger.info('Applying new conformation to clone')
        new_conf.save(conf_clone.structure_file) # overwrite the clone's structure with the new conformer
        
def perform_prop_analysis(polymer : Polymer, poly_logger : Logger, sim_dir_filters : Iterable[SimDirFilter], traj_sample_interval : int=1) -> None:
    '''Analyze trajectories to obtain polymer property data in reusable CSV form'''
    for sim_dir, (sim_paths, sim_params) in polymer.filter_sim_dirs(conditions=sim_dir_filters).items():
        poly_logger.info(f'Found trajectory {sim_paths.trajectory}')
        traj = polymer.load_traj(sim_dir)

        # save and plot RDF data
        poly_logger.info('Calculating pairwise radial distribution functions')
        rdf_dataframe = trajectory.acquire_rdfs(traj, max_rad=1.0*nanometer)
        rdf_save_path = sim_dir/'rdfs.csv'
        sim_paths.spatial_data = rdf_save_path
        rdf_dataframe.to_csv(rdf_save_path, index=False)

        # save and plot property data
        poly_logger.info('Calculating polymer shape properties')
        prop_dataframe = trajectory.acquire_time_props(traj, time_points=sim_params.time_points[::traj_sample_interval]) 
        prop_save_path = sim_dir/'time_series.csv'
        sim_paths.time_data = prop_save_path
        prop_dataframe.to_csv(prop_save_path, index=False)

        sim_paths.to_file(polymer.simulation_paths[sim_dir]) # update references to analyzed data files in path file
        poly_logger.info(f'Successfully exported trajectory analysis data')