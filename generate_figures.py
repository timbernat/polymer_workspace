'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
import re
from pathlib import Path
from collections import defaultdict
from itertools import combinations

# Numeric Imports
from math import ceil
import numpy as np
import pandas as pd
from openmm.unit import nanosecond

# Logging
import logging
logging.basicConfig(level=logging.INFO)
import matplotlib.pyplot as plt

from polysaccharide import LOGGERS_MASTER
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
from polysaccharide.graphics import plotutils
from polysaccharide.representation import Polymer, PolymerManager, has_sims
from polysaccharide.analysis import trajectory, statistics, equilibrium

# Static Paths
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
COLL_PATH = Path('Collections')

RESOURCE_PATH = Path('resources')
CHG_PARAM_PATH = RESOURCE_PATH / 'chg_templates'
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

src_coll_path = COLL_PATH / 'water_soluble_large'
save_dir = Path('figures_for_paper/100_ns_trials')
compare_dir = Path('colina_data')

is_long_sim = lambda sim_paths, sim_params : (sim_params.total_time == 100*nanosecond)
has_binary_traj = lambda sim_paths, sim_params : (sim_params.report_to_pdb == False)

CHG_CVTRS = {
    'Charge heatmap 2D' : 'SMARTS',
    'Charge heatmap 3D' : 'CXSMARTS'
}
equil_det = equilibrium.BinSegEquilDetector()

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    # load manager and create save folder and sets of structure converters and labels
    mgr = PolymerManager(src_coll_path)
    save_dir.mkdir(exist_ok=True)

    @mgr.logging_wrapper(loggers=loggers, proc_name='Figure generation', filters=has_sims)
    def generate_figures(polymer : Polymer, main_logger : logging.Logger, compare_dir : Path, obs_precision : int=3):
        '''Iterate through simulated polymers and generate data and plots'''
        # create folders to save to
        mol_name = polymer.mol_name
        parent_dir = save_dir / mol_name
        parent_dir.mkdir(exist_ok=True)
        observ = defaultdict(defaultdict)

        fig_time, ax_time = None, None
        fig_acf, ax_acf = None, None

        # Plot relative charges
        main_logger.info('Plotting charge-comparison heatmaps')
        chg_methods = polymer.charges.keys()
        for chg_pair in combinations(chg_methods, 2): # Plot all possible combinations of pairwise charge comparisons
            for label, cvtr_name in CHG_CVTRS.items():
                main_logger.info(f'Plotting {label} heatmap (via the "{cvtr_name}" converter)')
                fig_chg, ax_chg = polymer.compare_charges(*sorted(chg_pair), converter=cvtr_name)
                fig_chg.savefig(parent_dir / f'{label} - {mol_name}', bbox_inches='tight')
                plt.close()
        main_logger.info('') # palate cleanser
            
        # Plot observables
        for sim_dir, (sim_paths, sim_params) in polymer.filter_sim_dirs(conditions=(is_long_sim, has_binary_traj)).items():
            # read simulation trajectory data
            chg_method = sim_params.charge_method
            time_data = pd.read_csv(sim_paths.time_data)
            spatial_data = pd.read_csv(sim_paths.spatial_data)
            main_logger.info(f'Found trajectory data based on "{chg_method}" charges')

            # Plot RDFs
            main_logger.info('Plotting radial distribution functions (RDFs)')
            radii, rdfs = trajectory.rdfs_to_plot_data(spatial_data)
            num_rdfs = len(rdfs.columns)
            nrows_rdf = num_rdfs // 3
            ncols_rdf = ceil(num_rdfs / nrows_rdf)
            fig_rdf, ax_rdf = plotutils.plot_df_props(radii, rdfs, nrows=nrows_rdf, ncols=ncols_rdf)
            main_logger.info('Saving RDF plots')
            fig_rdf.savefig(parent_dir / f'RDFs - {mol_name}.png', bbox_inches='tight')
            plt.close()

            # Plot time series' and autocorrelation functions
            main_logger.info(f'Acquiring polymer property time series from {chg_method}-charged trajectory')
            times, time_series = trajectory.props_to_plot_data(time_data)
            if ax_time is None: # create new axes for first time...
                fig_time, ax_time = plotutils.plot_df_props(times, time_series)
            else: # otherwise, overlay plots
                for ax, (prop_name, data) in zip(ax_time, time_series.items()):
                    ax.plot(times, data)
                    ax.legend(chg_methods)

            main_logger.info(f'Acquiring polymer property autocorrelation functions (ACFs) from {chg_method}-charged trajectory')
            time_series_acfs = time_series.apply(statistics.autocorrelate)
            time_series_acfs = time_series_acfs.rename(mapper=lambda name : re.sub('\(.*?\)', 'ACF', name), axis='columns')
            if ax_acf is None: # create new axes for first time...
                fig_acf, ax_acf = plotutils.plot_df_props(times, time_series_acfs)
            else: # otherwise, overlay plots
                for ax, (prop_name, data) in zip(ax_acf, time_series_acfs.items()):
                    ax.plot(times, data)
                    ax.legend(chg_methods)

            # Determine equilibration times and average observables
            main_logger.info(f'Determining equilibration times for charge method {chg_method} trajectory')
            nrows, ncols = 1, len(time_series.columns)
            fig_equil, ax_equil = plotutils.presize_subplots(nrows, ncols)
            for ax, (prop_name, series) in zip(ax_equil, time_series.items()):
                series = np.array(series)
                equil_idx = equil_det.equil_loc(series)
                equil_time = times.to_numpy().flatten()[equil_idx]
                main_logger.info(f'Determined "{prop_name}" to have equilibrated after time {equil_time}') # TODO : find way to incorporate units into this

                obs = round(np.mean(series[equil_idx:]), obs_precision)
                observ[f'OpenFF 2.0.0 - {chg_method}'][prop_name] = obs
                main_logger.info(f'Determined equilibrium average of "{prop_name}" to be {obs}') 

                main_logger.info(f'Plotting equilibration cutoff curve for "{prop_name}"')
                ax.plot(times, series)
                ax.axvline(equil_time, color='r')
                ax.set_xlabel(times.columns[0])
                ax.set_ylabel(prop_name)
                ax.legend([f'<O> = {obs} (after t={equil_time})'], handlelength=0, loc='best') # position annotation in empty region
            main_logger.info('Saving equilibration plots')
            fig_equil.savefig(parent_dir / f'Equil times {chg_method} - {mol_name}.png', bbox_inches='tight')
            plt.close()
            main_logger.info('') # palate cleanser

        # Plot overlaid property and ACF data
        main_logger.info(f'Saving polymer property time series plots for {mol_name}')
        fig_time.legend()
        fig_time.savefig(parent_dir / f'Polymer Properties - {mol_name}.png', bbox_inches='tight')
        plt.close()
        
        main_logger.info(f'Saving polymer property ACF plots for {mol_name}')
        fig_acf.legend()
        fig_acf.savefig(parent_dir / f'Property ACFs - {mol_name}.png', bbox_inches='tight')
        plt.close()

        # Save CSV tables of observables, with experimental data appended
        main_logger.info(f'Acquiring experiemntal / literature observables for comparison')
        compare_path = compare_dir / f'{mol_name}_exp.csv'
        observ_df = pd.DataFrame(observ)
        compare_df = pd.read_csv(compare_path, index_col=0) # gets rid of unnecessary index col, enables smooth concatenation

        main_logger.info(f'Saving table of observables')
        observ_df_aug = pd.concat([observ_df, compare_df], axis=1)
        observ_df_aug.to_csv(parent_dir/ f'Observables - {mol_name}.csv')

    generate_figures(main_logger, compare_dir=compare_dir)