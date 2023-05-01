'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
from pathlib import Path
from openmm.unit import nanometer
from shutil import copyfile, rmtree

# Logging
import logging
logging.basicConfig(level=logging.INFO)

from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler

main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.solvation.solvents import WATER_TIP3P
from polysaccharide.molutils.polymer.building import build_linear_polymer_limited
from polysaccharide.charging.application import CHARGER_REGISTRY, ChargingParameters

# Static Paths
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
COLL_PATH = Path('Collections')

RESOURCE_PATH = Path('resources')
CHG_PARAM_PATH = RESOURCE_PATH / 'chg_templates'
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

def obtain_partial_charges(polymer : Polymer, main_logger : logging.Logger, chg_params : ChargingParameters) -> None:
    '''Ensure a Polymer has all partial charge sets'''
    # 0) LOAD MOLECULE AND TOPOLOGY, ATTEMPT TO APPLY LIBRARY CHARGES
    if not polymer.has_monomer_data:
        raise FileExistsError(f'No monomer JSONs found for {polymer.mol_name}')

    # 1) ENSURING CHARGES AND RELATED FILES FOR ALL CHARGING METHODS EXIST
    for chg_method in chg_params.charge_methods:
        chgr = CHARGER_REGISTRY[chg_method]()
        if chg_method == 'ABE10_averaged': # !NOTE! - critical that this not be the first key in the registry (has nothing to average over from scratch)
            residue_charges = polymer.residue_charges(
                averaging_charge_method=chg_params.averaging_charge_method,
                overwrite_charged_monomer_file=chg_params.overwrite_chg_mono
            )
            chgr.set_residue_charges(residue_charges)
        polymer.assert_charges_for(chgr, return_cmol=False)

    if (polymer.ff_file is None) or chg_params.overwrite_ff_xml: # can only reach if a charged monomer json already exists
        main_logger.info('Acquiring Force Field file with Library Charges')
        forcefield, lib_chgs = polymer.create_FF_file(xml_src=chg_params.base_ff_path, return_lib_chgs=True)

# ------------------------------------------------------------------------------

# Set parameters here
# source structures and reduced chain parameters
clear_prior = True 
source_coll_name = 'water_soluble_polymers'
chain_lim = 180
flip_term_group_labels = ['paam_modified']

# collection and charging paramters for reduced WaSPs
reduced_coll_name = 'water_soluble_reduced'
reduced_chg_params_path = CHG_PARAM_PATH / 'standard_chg_params.json'

# Collection, solvent, and simulation parameters for enlarged WaSPs
large_coll_name = 'water_soluble_large'
large_chg_params_path = CHG_PARAM_PATH / 'long_chain_chg_params.json'

solv_template    = RESOURCE_PATH/'inp_templates'/'solv_polymer_template_box.inp'
desired_solvent = WATER_TIP3P
exclusion = 1*nanometer

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main_logger.info('STAGE (0) : Creating paths and empty directories')
    # source polymers
    poly_source_path = COMPAT_PDB_PATH / source_coll_name
    structure_path   = poly_source_path / f'{poly_source_path.name}_structures'
    monomer_path     = poly_source_path / f'{poly_source_path.name}_monomers'

    src_coll_path = COLL_PATH / source_coll_name
    # don;t need to mkdir these paths (should already exist)

    # reduced polymers
    reduced_src_path   = COMPAT_PDB_PATH / reduced_coll_name
    reduced_structure_path = reduced_src_path / f'{reduced_src_path.name}_structures'
    reduced_monomer_path   = reduced_src_path / f'{reduced_src_path.name}_monomers'

    reduced_coll_path = COLL_PATH / reduced_coll_name

    # full-size polymers
    large_coll_path = COLL_PATH / large_coll_name

    # optional clean-up
    if clear_prior:
        main_logger.warning('Cleaning up directories from previous sessions')
        for dir in (reduced_src_path, src_coll_path, reduced_coll_path, large_coll_path):
            if dir.is_dir() and dir.exists():
                rmtree(dir)

    reduced_src_path.mkdir(  exist_ok=True)
    reduced_monomer_path.mkdir(  exist_ok=True)
    reduced_structure_path.mkdir(exist_ok=True)
    large_coll_path.mkdir(exist_ok=True)

# STRUCTURE GENERATION
    main_logger.info('STAGE (1) : Generating reduced-size structures and monomers for WaSPs')

    mgr_src = PolymerManager(src_coll_path)
    if not mgr_src.polymers: # ensure originals have been created
        mgr_src.populate_collection(struct_dir=structure_path, monomer_dir=monomer_path) 
    
    # generate new structures
    reverse = False # needed
    for pdir in mgr_src.polymers_list:
        if pdir.solvent is None: # explicitly don't copy solvated systems
            print(pdir.mol_name)
            monomers = pdir.monomer_data['monomers']
            if pdir.mol_name in flip_term_group_labels:
                monomers.pop('paam_SPECIAL_TERM')   # special case needed for extra messy term group in PAAM...
                reverse = True                      # ...TODO : find generalized way to isolate head and tail groups from larger collection of monomers

            chain = build_linear_polymer_limited(monomers, max_chain_len=chain_lim, reverse_term_labels=reverse)
            chain.save(str(reduced_structure_path/f'{pdir.mol_name}.pdb'), overwrite=True)
            copyfile(pdir.monomer_file, reduced_monomer_path/f'{pdir.mol_name}.json')
    
# COMPUTE ABE10 EXACT AND AVERAGED CHARGES FOR REDUCTIONS
    main_logger.info('STAGE (2) Populate and load collection from reduced structures, then charge and compute averaged residues')
    mgr_reduced = PolymerManager(reduced_coll_path)

    if not mgr_reduced.polymers: # will be empty if not yet instantiated or if reset prior
        mgr_reduced.populate_collection(struct_dir=reduced_structure_path, monomer_dir=reduced_monomer_path) # don't solvate yet, as this is unnecessary

    reduced_chg_params = ChargingParameters.from_file(reduced_chg_params_path)
    with ProcessLogHandler(filedir=mgr_reduced.log_dir, loggers=loggers, proc_name='Charging of reductions', timestamp=True) as msf_handler:
        for i, (mol_name, polymer) in enumerate(mgr_reduced.polymers.items()):
            main_logger.info(f'Current molecule: "{mol_name}" ({i + 1}/{mgr_reduced.n_mols})') # +1 converts to more human-readable 1-index for step count
            with msf_handler.subhandler(filedir=polymer.logs, loggers=loggers, proc_name='Charging', timestamp=True) as subhandler: # also log actions to individual Polymers
                obtain_partial_charges(polymer, main_logger, reduced_chg_params)

# TRANSFER AVERAGED CHARGES TO FULL-SIZE STRUCTURES
    main_logger.info('STAGE (3) Transfer averaged charges from reduced WaSPs to full-sized structures ')
    mgr_large = PolymerManager(large_coll_path)

    with ProcessLogHandler(filedir=mgr_reduced.log_dir, loggers=loggers, proc_name='Charge transfer'):
        for polymer in mgr_reduced.polymers_list:
            if polymer.charges and polymer.has_monomer_data_charged:
                polymer.clone(
                    dest_dir=large_coll_path,
                    clone_name=polymer.base_mol_name, # exclude solvent (will need to resolvate with new structure later)
                    clone_solvent=False,
                    clone_structures=False,
                    clone_monomers=True, # keep only charged monomer information
                    clone_ff=False,
                    clone_charges=False,
                    clone_sims= False
                )

    # copy large structures over to clones, initialize as PolymerManager
    with ProcessLogHandler(filedir=mgr_large.log_dir, loggers=loggers, proc_name='Large structure alignment'):
        mgr_large.update_collection() # update collection to reflect newly created enlarged copies
        for polymer in mgr_large.polymers_list: # copy large structures over
            polymer.populate_pdb(structure_path)
            polymer.solvate(solvent=desired_solvent, template_path=solv_template, exclusion=exclusion) # assumes only 1 solvent
            
        for polymer in mgr_large.polymers_list:
            if polymer.solvent == desired_solvent:
                polymer.create_charged_monomer_file(residue_charges=polymer.monomer_data_charged['charges']) # needed to incorporate solvent charges into charged JSON
    
    print(mgr_large.polymers.keys())

    large_chg_params = ChargingParameters.from_file(large_chg_params_path)
    with ProcessLogHandler(filedir=mgr_large.log_dir, loggers=loggers, proc_name='Charging of large chains', timestamp=True) as msf_handler:
        for i, polymer in enumerate(mgr_large.polymers_list):
            main_logger.info(f'Current molecule: "{polymer.mol_name}" ({i + 1}/{mgr_large.n_mols})') # +1 converts to more human-readable 1-index for step count
            with msf_handler.subhandler(filedir=polymer.logs, loggers=loggers, proc_name='Charging', timestamp=True) as subhandler: # also log actions to individual Polymers
                obtain_partial_charges(polymer, main_logger, large_chg_params)
    
# SIMULATE AND COMPUTE PROPERTIES
    # main_logger.info('STAGE (4) Run simulation and analysis over large WaSPs, obtain final trajectories and data')