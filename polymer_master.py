'''Master script for accessing and deploying polymer workflow components'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
from pathlib import Path

# Polymer Imports
from polysaccharide.polymer.management import PolymerManager, MolFilterBuffer
import polysaccharide.workflow.components as components


# CLI argument parsing
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-src', '--source_path' , help='The Path to the target collection of Polymers', required=True, type=Path)
parser.add_argument('-ps', '--parallelize_sbatch', help='Whether or not to dispatch jobs in parallel via a job submission script', action='store_true')
parser.add_argument('-sb', '--sbatch_script'     , help='Name of the target slurm job script to use for submission', default='slurm_dispatch.job', type=Path)
parser.add_argument('-jid', '--collect_job_ids'  , help='Whether or not to gather job IDs when submitting (useful for creating dependencies in serial workflows)', action='store_true')
subparsers = parser.add_subparsers(help='Concrete implementations of various workflow components', dest='component')

for comp_name, Component in components.WorkflowComponent.registry.items():
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
    Component = components.WorkflowComponent.registry[args.component]

    # Defining filters and component protocols
    comp = Component.from_argparse(args)
    filters = comp.assert_filter_prefs(molbuf)

    # branch for parallelization
    if args.parallelize_sbatch:
        slurm_comp = components._SlurmSbatch(
            comp,
            sbatch_script=args.sbatch_script,
            python_script_name=__file__, # recall this script for singleton molecule
            source_path=args.source_path,
            collect_job_ids=args.collect_job_ids
        )

        task_fn = src_mgr.logging_wrapper(proc_name=f'{Component.__name__}Dispatch', filters=filters)(slurm_comp.make_polymer_fn()) 
    else:
        task_fn = src_mgr.logging_wrapper(proc_name=Component.__name__, filters=filters)(comp.make_polymer_fn()) 
        
    # Execution
    task_fn()
    if args.parallelize_sbatch and args.collect_job_ids:
        print(slurm_comp.dependency_str)
    