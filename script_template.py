'''Generic template for defining CLI-executable polymer action scripts'''

'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

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

# Polymer Imports
from polysaccharide.representation import Polymer, PolymerManager
from polysaccharide.representation import identity, filter_factory_by_attr
from polysaccharide.charging.application import ChargingParameters

# Static Paths
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
RESOURCE_PATH = impres.files(resources)

# CLI arg parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=''
)
parser.add_argument('-src', '--source_name'  , help='The name of the target collection of Polymers', required=True)

args = parser.parse_args()

# Arg processing
# ------------------------------------------------------------------------------

source_path = COLL_PATH / args.source_name
...
filters = [identity]

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    mgr = PolymerManager(source_path)

    @mgr.logging_wrapper(loggers, proc_name='...', filters=filters)
    def execute(polymer : Polymer, *args, **kwargs):
        pass

    execute()