'''Master script for accessing and deploying polymer workflow components'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path

# Polymer Imports
from polysaccharide.polymer.management import PolymerManager, MolFilterBuffer
import polymer_workflow


# CLI argument parsing
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-src', '--source_path' , help='The Path to the target collection of Polymers', required=True, type=Path)
subparsers = parser.add_subparsers(help='Concrete implementations of various workflow components', dest='component')

for comp_name, Component in polymer_workflow.WorkflowComponent.registry.items():
    subparser = subparsers.add_parser(comp_name, help=Component.desc)
    Component.argparse_inject(subparser) # generate dedicated subparser for each component
    MolFilterBuffer.argparse_inject(subparser) # add filtering options 

args = parser.parse_args()

# Compilation and execution
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    # Processing args into principal objects
    src_mgr = PolymerManager(args.source_path)
    molbuf = MolFilterBuffer.from_argparse(args)
    Component = polymer_workflow.WorkflowComponent.registry[args.component]

    # Defining filters and component protocols
    comp = Component.from_argparse(args)
    filters = comp.assert_filter_prefs(molbuf)
    print(comp.__dict__, filters)

    task_fn = src_mgr.logging_wrapper(proc_name=Component.__name__, filters=filters)(comp.make_polymer_fn()) 
    
    # Execution
    task_fn()
