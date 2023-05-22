'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Logging
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO, force=True)

from polysaccharide import LOGGERS_MASTER
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
import importlib_resources as impres
import resources

from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import is_unsolvated, is_uncharged, filter_factory_by_attr
from polysaccharide.charging.application import ChargingParameters

# Static Paths
COLL_PATH = Path('Collections')
RESOURCE_PATH = Path('resources')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

src_coll_name = 'simple_polymers'
chg_params_name = 'standard_chg_params.json'

# src_coll_name = 'water_soluble_large_conf'
# chg_params_name = 'long_chain_chg_params.json'

targ_mols = ('naturalrubber', 'polyvinylchloride')
desired_mol = filter_factory_by_attr('mol_name', condition=lambda name : name in targ_mols)

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    src_coll_path = COLL_PATH / src_coll_name
    mgr = PolymerManager(src_coll_path)
    
    generate_charge_files = mgr.logging_wrapper(
        loggers,
        proc_name=f'Charge assignment {src_coll_name}',
        filters=(is_uncharged, desired_mol) #, is_unsolvated)
    )(Polymer.obtain_partial_charges) 

    with impres.path(resources.chg_templates, chg_params_name) as chg_params_path:
        chg_params = ChargingParameters.from_file(chg_params_path)
        generate_charge_files(chg_params)
