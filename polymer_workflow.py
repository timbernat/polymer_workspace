# Logging
import logging

# Typing and Subclassing
from typing import Any, Callable, ClassVar, Iterable, Optional, Union
from abc import ABC, abstractmethod, abstractproperty

# Generic imports
from pathlib import Path
from shutil import copyfile
from time import sleep
from abc import abstractstaticmethod

from argparse import ArgumentParser, Namespace

# Resource imports
import importlib_resources as impres
import resources

avail_chg_templates = ', '.join(
    path.name
        for path in resources.AVAIL_RESOURCES['chg_templates']
)

avail_sim_templates = ', '.join(
    path.name
        for path in resources.AVAIL_RESOURCES['sim_templates']
)

# Polymer Imports
from polysaccharide.solvation import solvents as all_solvents

from polysaccharide.charging.application import ChargingParameters, CHARGER_REGISTRY
from polysaccharide.charging.averaging import get_averaged_residue_charges, AveragingCharger

from polysaccharide.simulation.records import SimulationParameters
from polysaccharide.simulation.execution import run_simulation

from polysaccharide.polymer.representation import Polymer
from polysaccharide.polymer.management import PolymerManager, PolymerFunction, MolFilterBuffer
from polysaccharide.polymer.filtering import MolFilter, has_sims, has_monomers_chgd, is_base
from polysaccharide.polymer.filtering import SimDirFilter, has_binary_traj

from polysaccharide.polymer.monomer import estimate_max_DOP, estimate_chain_len
from polysaccharide.polymer.building import build_linear_polymer
from polysaccharide.polymer.exceptions import ExcessiveChainLengthError

from polysaccharide.analysis import trajectory

# Molecular Dynamics
from openmm.unit import nanosecond # time
from openmm.unit import nanometer, angstrom # length


# Base class
class WorkflowComponent(ABC):
    @abstractproperty
    @classmethod
    def desc(self) -> str:
        '''Brief description to accompany component'''
        ...
    
    @abstractproperty
    @classmethod
    def name(self) -> str:
        '''Brief name to label component'''
        ...

    @abstractmethod
    def __init__(self, *args, **kwargs):
        ...

    @abstractstaticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        '''Flexible support for instantiating addition to argparse in an existing script'''
        ...

    @abstractmethod
    def make_polymer_fn(self) -> PolymerFunction:
        ...

    @classmethod
    def from_argparse(cls, args : Namespace) -> 'WorkflowComponent':
        '''Initialize from an argparse Namespace'''
        return cls(**vars(args))

    def assert_filter_prefs(self, molbuf : MolFilterBuffer) -> list[MolFilter]:
        '''Assert any additional preferences for filters beyond the default molecule filters'''
        return molbuf.filters # default to base filters

    @classmethod
    @property
    def registry(cls) -> dict[str, 'WorkflowComponent']:
        '''Name-indexed dict of all inherited Component implementations'''
        return {
            subcomp.name : subcomp
                for subcomp in cls.__subclasses__()
        }


# Concrete child class implementations
class DummyCalculation(WorkflowComponent):
    desc = 'Computes RDF and property time series data and saving to csvs for plotting and analysis'
    name = 'dummy'

    def __init__(self, wait_time : int, **kwargs):
        '''Initialize wait time to simulate non-trivial task'''
        self.wait_time = wait_time

    @staticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        '''Flexible support for instantiating addition to argparse in an existing script'''
        parser.add_argument('-w', '--wait_time', help='Number of seconds for dummy task to wait', action='store', type=int)

    def make_polymer_fn(self) -> PolymerFunction:
        '''Create wrapper for handling in logger'''
        def polymer_fn(polymer : Polymer, poly_logger : logging.Logger) -> None:
            '''Dummy function for testing script dispatch'''
            poly_logger.info(f'Performing fake calculation for {polymer.mol_name} for {self.wait_time} seconds')
            sleep(self.wait_time)

        return polymer_fn

class ChargeAssignment(WorkflowComponent):
    desc = 'Partial charge assignment'
    name = 'charge'

    def __init__(self, chg_params_name : str, **kwargs):
        '''Load charging parameters from central resource file'''
        chg_params_path = impres.files(resources.chg_templates) / chg_params_name
        if not chg_params_path.suffix:
            chg_params_path = chg_params_path.with_name(f'{chg_params_path.stem}.json') # ensure charge params path has correct extension

        self.chg_params = ChargingParameters.from_file(chg_params_path)

    @staticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        parser.add_argument('-cp', '--chg_params_name', help=f'Name of the charging parameters preset file to load for charging (available files are {avail_chg_templates})', required=True)

    # TOSELF : overwrite / charge status force in assert_filter_prefs?

    def make_polymer_fn(self) -> PolymerFunction:
        '''Create wrapper for handling in logger'''
        def polymer_fn(polymer : Polymer, poly_logger : logging.Logger) -> None:
            '''Ensure a Polymer has all partial charge sets'''
            assert(polymer.has_monomer_info)

            # 1) ENSURING CHARGES AND RELATED FILES FOR ALL BASE CHARGING METHODS EXIST
            for chg_method in self.chg_params.charge_methods:
                chgr = CHARGER_REGISTRY[chg_method]()
                polymer.assert_charges_for(chgr, return_cmol=False)
                poly_logger.info('') # add gaps between charge method for breathing room
                
            # 2) GENERATE RESIDUE-AVERAGED LIBRARY CHARGES FROM THE CHARGE SET OF CHOICE
            if polymer.has_monomer_info_charged: # load precomputed residue charges if existing set is found
                poly_logger.info(f'Found charged JSON, loading averaged charges for {polymer.mol_name} residues from file')
                residue_charges = polymer.monomer_info_charged.charges
            else: # otherwise, recompute them from the charge set of choice
                poly_logger.info(f'Averaging {self.chg_params.averaging_charge_method} charges over {polymer.mol_name} residues')
                residue_charges = get_averaged_residue_charges(
                    cmol=polymer.charged_offmol(self.chg_params.averaging_charge_method), # TODO : add check to ensure this ISN'T and averaging method
                    monomer_info=polymer.monomer_info_uncharged
                )
            
            if (not polymer.has_monomer_info_charged) or self.chg_params.overwrite_chg_mono: # cast charges residues to file if none exist or is explicitly overwriting
                polymer.create_charged_monomer_file(residue_charges) # TOSELF : this is separate from the above clauses as it might be called regardless of existing charges during overwrite

            avg_chgr = AveragingCharger() # generate precomputed charge set for full molecule from residue library charges
            avg_chgr.set_residue_charges(residue_charges)
            polymer.assert_charges_for(avg_chgr, return_cmol=False)
        
        return polymer_fn
    
class TrajectoryAnalysis(WorkflowComponent):
    desc = 'Computes RDF and property time series data and saving to csvs for plotting and analysis'
    name = 'analyze'

    def __init__(self, sim_time : Optional[int], traj_sample_interval : int=1, **kwargs):
        '''Defining simulation-based filters'''
        self.sim_dir_filters = [has_binary_traj]
        if sim_time is not None:
            is_long_sim = lambda sim_paths, sim_params : (sim_params.total_time == sim_time*nanosecond)
            self.sim_dir_filters.append(is_long_sim)
        
        self.traj_sample_interval = traj_sample_interval

    @staticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        '''Flexible support for instantiating addition to argparse in an existing script'''
        parser.add_argument('-t', '--sim_time'              , help='If set, will only analyze trajectories run for this number of nanoseconds', action='store', type=float)
        parser.add_argument('-tsi', '--traj_sample_interval', help='How often to sample trajectory frames when loading (equilvalent to "stride" in mdtraj); useful for huge trajectories', action='store', type=int, default=1)

    def assert_filter_prefs(self, molbuf : MolFilterBuffer) -> list[MolFilter]:
        '''Assert any additional preferences for filters beyond the default molecule filters'''
        mol_filters = molbuf.filters
        mol_filters.append(has_sims)

        return mol_filters

    def make_polymer_fn(self) -> PolymerFunction:
        '''Create wrapper for handling in logger'''
        def polymer_fn(polymer : Polymer, poly_logger : logging.Logger) -> None:
            '''Analyze trajectories to obtain polymer property data in reusable CSV form'''
            for sim_dir, (sim_paths, sim_params) in polymer.filter_sim_dirs(conditions=self.sim_dir_filters).items():
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
                prop_dataframe = trajectory.acquire_time_props(traj, time_points=sim_params.time_points[::self.traj_sample_interval]) 
                prop_save_path = sim_dir/'time_series.csv'
                sim_paths.time_data = prop_save_path
                prop_dataframe.to_csv(prop_save_path, index=False)

                sim_paths.to_file(polymer.simulation_paths[sim_dir]) # update references to analyzed data files in path file
                poly_logger.info(f'Successfully exported trajectory analysis data')
            
        return polymer_fn
    
class BuildReducedStructures(WorkflowComponent):
    desc = 'Creates directory of reduced-chain-length structure PDBs and monomer JSONs from a collection of linear Polymers'
    name = 'redux'

    def __init__(self, struct_output : Path, mono_output : Path, DOP : Optional[int]=None, max_chain_len : Optional[int]=None, chain_len_limit : int=300, flip_term_labels : Optional[Iterable]=None):
        '''Initialize sizes of new chains to build, along with locations to output structure files to'''
        self.struct_output = struct_output  
        self.mono_output = mono_output 

        if not (max_chain_len or DOP):
            raise ValueError('Must provide EITHER a maximum chain length OR a degree of polymerization (provided neither)')

        if max_chain_len and DOP:
            raise ValueError('Must provide EITHER a maximum chain length OR a degree of polymerization (provided both)')
        
        self.max_chain_len = max_chain_len  
        self.DOP = DOP  
        self.chain_len_limit = chain_len_limit  

        self.flip_term_labels = flip_term_labels  

    @staticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        '''Flexible support for instantiating addition to argparse in an existing script'''
        parser.add_argument('-pdb', '--struct_output'  , help='The name of the directory to output generated PDB structure to', type=Path)
        parser.add_argument('-mono' , '--mono_output'  , help='The name of the directory to output generated JSON monomer files to', type=Path)
        parser.add_argument('-N', '--max_chain_len'    , help='Maximum number of atoms in any of the reduced chain generated. If this is specified, CANNOT specify DOP', type=int)
        parser.add_argument('-D', '--DOP'              , help='The number of monomer units to include in the generated reductions.  If this is specified, CANNOT specify max_chain_len', type=int)
        parser.add_argument('-lim', '--chain_len_limit', help='The maximum allowable size for a chain to be built to; any chains attempted to be built larger than this limit will raise an error', type=int, default=300)
        parser.add_argument('-f', '--flip_term_labels' , help='Names of the chains on which to reverse the order of head/tail terminal group labels (only works for linear homopolymers!)', action='store', nargs='+', default=tuple())

    def assert_filter_prefs(self, molbuf : MolFilterBuffer) -> list[MolFilter]:
        '''Assert any additional preferences for filters beyond the default molecule filters'''
        molbuf.solvent = False # force preference for unsolvated molecules - makes logic for monomer selection cleaner (don;t need to worry about wrong number of monomers due to solvent)
        return molbuf.filters 

    def make_polymer_fn(self) -> PolymerFunction:
        '''Create wrapper for handling in logger'''
        def polymer_fn(polymer : Polymer, poly_logger : logging.Logger) -> None:
            '''Builds new PDB structures of the desired size from the monomers of an existing Polymer'''
            monomer_smarts = polymer.monomer_info.monomers # create copy to avoid popping from original

            if DOP: # NOTE : this only works as intended because of the exclusivity check during arg processing
                max_chain_len = estimate_chain_len(monomer_smarts, DOP)

            if max_chain_len:
                DOP = estimate_max_DOP(monomer_smarts, max_chain_len)
            
            if self.max_chain_len > self.chain_len_limit:
                raise ExcessiveChainLengthError(f'Cannot create reduction with over {self.chain_len_limit} atoms (requested {max_chain_len})')
            
            chain = build_linear_polymer(monomer_smarts, DOP=DOP, reverse_term_labels=(polymer.mol_name in self.flip_term_labels))
            chain.save(str(self.struct_output/f'{polymer.mol_name}.pdb'), overwrite=True)
            copyfile(polymer.monomer_file_uncharged, self.mono_output/f'{polymer.mol_name}.json')

        return polymer_fn
    
class RunSimulations(WorkflowComponent):
    desc = 'Prepares and integrates MD simulation for chosen molecules in OpenMM'
    name = 'simulate'

    def __init__(self, sim_param_names : Iterable[str], **kwargs):
        self.all_sim_params= []
        for sim_param_name in sim_param_names:
            sim_param_path = impres.files(resources.sim_templates)/ sim_param_name
            if not sim_param_path.suffix:
                sim_param_path = sim_param_path.with_name(f'{sim_param_path.stem}.json') # ensure charge params path has correct extension

            self.all_sim_params.append( SimulationParameters.from_file(sim_param_path) )

    @staticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        '''Flexible support for instantiating addition to argparse in an existing script'''
        parser.add_argument('-sp', '--sim_param_names', help=f'Name of the simulation parameters preset file(s) to load for simulation (available files are {avail_sim_templates})', action='store', nargs='+', required=True)

    def assert_filter_prefs(self, molbuf : MolFilterBuffer) -> list[MolFilter]:
        '''Assert any additional preferences for filters beyond the default molecule filters'''
        molbuf.charges = True # force preference for charged molecules (can't run simulations otherwise)
        return molbuf.filters

    def make_polymer_fn(self) -> PolymerFunction:
        '''Create wrapper for handling in logger'''
        def polymer_fn(polymer : Polymer, poly_logger : logging.Logger) -> None:
            '''Run OpenMM simulation(s) according to sets of predefined simulation parameters'''
            N = len(self.all_sim_params)
            for i, sim_params in enumerate(self.all_sim_params):
                poly_logger.info(f'Running simulation {i + 1} / {N}')
                interchange = polymer.interchange(
                    forcefield_path=sim_params.forcefield_path,
                    charge_method=sim_params.charge_method,
                    periodic=sim_params.periodic
                )

                sim_folder = polymer.make_sim_dir()
                run_simulation(interchange, sim_params=sim_params, output_folder=sim_folder, output_name=polymer.mol_name)

        return polymer_fn
    
class Solvate(WorkflowComponent):
    desc = 'Solvate molecules in sets of 1 or more desired solvents'
    name = 'solvate'

    def __init__(self, solvents : Iterable[str], template_name : Path, exclusion : float=None, **kwargs):
        if not solvents:
            raise ValueError('Must specify at least 1 solvent')
        self.solvents = [
            getattr(all_solvents, solvent_name)
                for solvent_name in solvents
        ]
        self.template_path = impres.files(resources.inp_templates) / template_name
        self.exclusion = exclusion * nanometer # assign units

    @staticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        '''Flexible support for instantiating addition to argparse in an existing script'''
        parser.add_argument('-s', '--solvents'     , help='Names of all solvent molecule to solvate the target systems in' , action='store', nargs='+', default=['WATER_TIP3P'])
        parser.add_argument('-t', '--template_name', help='Name of the packmol input template file to use for solvation', action='store', default='solv_polymer_template_box.inp')
        parser.add_argument('-e', '--exclusion'    , help='Distance (in nm) between the bounding box of the molecule and the simiulation / solvation box', action='store', type=float, default=1.0)

    def assert_filter_prefs(self, molbuf : MolFilterBuffer) -> list[MolFilter]:
        '''Assert any additional preferences for filters beyond the default molecule filters'''
        molbuf.solvent = False # force preference for only unsolvated molecules (don't want to attempt solvation twice)
        return molbuf.filters

    def make_polymer_fn(self) -> PolymerFunction:
        '''Create wrapper for handling in logger'''
        def polymer_fn(polymer : Polymer, poly_logger : logging.Logger) -> None:
            '''Fill a box around a polymer with solvent'''
            polymer.solvate(self.solvents, template_path=self.template_path, exclusion=self.exclusion)
 
        return polymer_fn
    
class TransferMonomerCharges(WorkflowComponent):
    desc = 'Transfer residue-averaged charges from reductions to full-sized collections'
    name = 'transfer'

    def __init__(self, target_path : str, **kwargs):
        self.charger_mgr = PolymerManager(target_path)

    @staticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        '''Flexible support for instantiating addition to argparse in an existing script'''
        parser.add_argument('-targ', '--target_path', help='The path to the target output collection of Polymers to move charged monomers to', type=Path, required=True)

    def assert_filter_prefs(self, molbuf : MolFilterBuffer) -> list[MolFilter]:
        '''Assert any additional preferences for filters beyond the default molecule filters'''
        molbuf.charges = True # force preference for charges 
        mol_filters = molbuf.filters
        mol_filters.append(has_monomers_chgd) # also assert that, not only do charges exist, but that they've been monomer-averaged

        return mol_filters

    def make_polymer_fn(self) -> PolymerFunction:
        '''Create wrapper for handling in logger'''
        def polymer_fn(polymer : Polymer, poly_logger : logging.Logger) -> None:
            '''Copies charged monomer files from a corresponding Polymer in another collection'''
            counterpart = self.charged_mgr.polymers[polymer.mol_name]
            assert(polymer.solvent == counterpart.solvent)

            if not polymer.has_monomer_info_uncharged:
                counterpart.transfer_file_attr('monomer_file_uncharged', polymer)
            
            if not polymer.has_monomer_info_charged:
                counterpart.transfer_file_attr('monomer_file_charged', polymer)
        
        return polymer_fn
    
class VacuumAnneal(WorkflowComponent): # TODO : decompose this into cloning, sim (already implemented), and structure transfer components
    desc = 'Generate arbitrary number varied starting structures via short, high-T MD simulations'
    name = 'anneal'

    def __init__(self, sim_params_name : str, num_new_confs : int, snapshot_idx : int=-1, **kwargs):
        '''Define parameters for vacuum anneal, along with number of new conformers'''
        sim_param_path = impres.files(resources.sim_templates) / sim_params_name
        if not sim_param_path.suffix:
            sim_param_path = sim_param_path.with_name(f'{sim_param_path.stem}.json') # ensure charge params path has correct extension
        self.sim_param_path = sim_param_path

        self.num_new_confs = num_new_confs
        self.snapshot_idx = snapshot_idx

    @staticmethod
    def argparse_inject(parser : ArgumentParser) -> None:
        '''Flexible support for instantiating addition to argparse in an existing script'''
        parser.add_argument('-sim', '--sim_params'  , help='Name of the simulation parameters preset file to load for simulation', default='vacuum_anneal.json')
        parser.add_argument('-r', '--num_replicates', help='Number of total conformers to generate (not counting the original)', type=int, default=4)
    
    def assert_filter_prefs(self, molbuf : MolFilterBuffer) -> list[MolFilter]:
        '''Assert any additional preferences for filters beyond the default molecule filters'''
        molbuf.solvent = False # force preference for unsolvated molecules (VACUUM anneal)
        molbuf.charges = True  # force preference for charged molecules (can;t run sim otherwise)
        mol_filters = molbuf.filters
        mol_filters.append(is_base) # also assert that only base molecules are annealed
        
        return mol_filters # default to base filters

    def make_polymer_fn(self) -> PolymerFunction:
        '''Create wrapper for handling in logger'''
        def polymer_fn(polymer : Polymer, poly_logger : logging.Logger) -> None:
            '''Run quick vacuum NVT sim at high T for specified number of runs to generate perturbed starting structures'''
            for i in range(self.num_new_confs):
                # generate clone to anneal
                conf_clone = polymer.clone(
                    clone_name=f'{polymer.base_mol_name}_conf_{i + 1}',
                    clone_solvent=True,
                    clone_structures=True,
                    clone_monomers=True,
                    clone_charges=True,
                    clone_sims=False
                )

                # runn simulation
                sim_step = RunSimulations(self.sim_param_path.stem)
                simulate_polymer = sim_step.make_polymer_fn()
                simulate_polymer(polymer, poly_logger)
                
                # replace clone's starting structure with new anneal structure
                poly_logger.info('Extracting final conformation from simulation')
                traj = polymer.load_traj(polymer.newest_sim_dir)
                new_conf = traj[self.snapshot_idx]
                poly_logger.info('Applying new conformation to clone')
                new_conf.save(conf_clone.structure_file) # overwrite the clone's structure with the new conformer
                
        return polymer_fn
    