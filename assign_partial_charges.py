'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Logging
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

from polysaccharide import LOGGERS_MASTER
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
import importlib_resources as impres
import resources

from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import is_unsolvated, is_charged
from polysaccharide.charging.application import CHARGER_REGISTRY, ChargingParameters

# Static Paths
COLL_PATH = Path('Collections')
RESOURCE_PATH = Path('resources')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

# src_coll_name = 'water_soluble_small_conf'
# chg_params_name = 'standard_chg_params.json'

src_coll_name = 'water_soluble_large_conf'
chg_params_name = 'long_chain_chg_params.json'

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    src_coll_path = COLL_PATH / src_coll_name
    mgr = PolymerManager(src_coll_path)

    @mgr.logging_wrapper(loggers, proc_name=f'Charge assignment {src_coll_name}') #, filters=(is_unsolvated))
    def obtain_partial_charges(polymer : Polymer, chg_params : ChargingParameters) -> None:
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
            forcefield = polymer.create_FF_file(xml_src=chg_params.base_ff_path, return_lib_chgs=True)

    with impres.path(resources.chg_templates, chg_params_name) as chg_params_path:
        chg_params = ChargingParameters.from_file(chg_params_path)
        obtain_partial_charges(chg_params)




 