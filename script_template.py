# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path

# Polymer Imports
from polysaccharide.polymer.management import PolymerManager, MolFilterBuffer
import polysaccharide.workflow.components as components

# CLI arg parsing
# ------------------------------------------------------------------------------

COMP_NAME = 'ChargeAssignment'
Component = getattr(components, COMP_NAME)

parser = argparse.ArgumentParser(description=Component.desc)
parser.add_argument('-src', '--source_path' , help='The Path to the target collection of Polymers', required=True, type=Path)
Component.argparse_inject(parser)
MolFilterBuffer.argparse_inject(parser)

# Arg processing
# ------------------------------------------------------------------------------

args = parser.parse_args()

src_mgr = PolymerManager(args.source_path)
comp = Component.from_argparse(args)
molbuf = MolFilterBuffer.from_argparse(args)

# Execution
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    task_fn = src_mgr.logging_wrapper(
        proc_name=Component.__name__,
        filters=molbuf.filters
    )(comp.make_polymer_fn()) 
    
    task_fn()
