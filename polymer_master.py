'''Master script for accessing and deploying polymer workflow components'''

# Logging
import logging
logging.basicConfig(level=logging.INFO, force=True)

# Generic imports
import argparse
import sys, subprocess

import re
from typing import Union
from pathlib import Path

def extract_job_id(shell_str : Union[str, bytes]) -> str:
    '''Extracts job id from shell-echoed string after slurm job submission'''
    JOB_ID_RE = re.compile(r'Submitted batch job (\d+)')
    if not isinstance(shell_str, str):
        shell_str = str(shell_str)

    return re.search(JOB_ID_RE, str(shell_str)).groups()[0]

# Polymer Imports
from polysaccharide.polymer.management import PolymerManager, MolFilterBuffer
import polymer_workflow


# CLI argument parsing
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-src', '--source_path' , help='The Path to the target collection of Polymers', required=True, type=Path)
parser.add_argument('-sb', '--sbatch_script', help='Whether or not to dispatch jobs in parallel via a job submission script')#, type=Path)
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

    if args.sbatch_script is None:
        # Execution
        task_fn = src_mgr.logging_wrapper(proc_name=Component.__name__, filters=filters)(comp.make_polymer_fn()) 
        task_fn()
    else:
        arg_start_idx = sys.argv.index(Component.name) + 1 # find where Component-specific arguments begin (immediately after job type spec)
        arg_str = ' '.join(sys.argv[arg_start_idx:])       # collate into passable string
        
        job_ids = []
        for (mol_name, polymer) in src_mgr.filtered_by(filters):
            cmd = ' '.join([
                'sbatch',
                f'--job-name "{mol_name}_dispatch"',
                f'--output "slurm_logs/{mol_name}_dispatch.log"',
                args.sbatch_script,
                str(args.source_path),
                Component.name,
                polymer.mol_name,
                arg_str
            ])
            print(cmd)
            # slurm_out = subprocess.check_output([cmd], shell=True)
            # job_ids.append(extract_job_id(slurm_out))
        print(job_ids)
    
