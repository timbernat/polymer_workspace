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

# Logging
import logging
logging.basicConfig(level=logging.INFO)
import matplotlib.pyplot as plt

from polysaccharide import LOGGERS_MASTER
from polysaccharide.logutils import ProcessLogHandler

main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER]

# Polymer Imports
from polysaccharide.graphics import plotutils
from polysaccharide.representation import PolymerManager, has_sims
from polysaccharide.simulation import SimulationParameters, SimulationPaths
from polysaccharide.analysis import trajectory, statistics

# Static Paths
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')
COLL_PATH = Path('Collections')

RESOURCE_PATH = Path('resources')
CHG_PARAM_PATH = RESOURCE_PATH / 'chg_templates'
SIM_PARAM_PATH = RESOURCE_PATH / 'sim_templates'

# ------------------------------------------------------------------------------

src_coll_path = COLL_PATH / 'water_soluble_large'
save_dir = Path('figures_for_paper')
compare_dir = Path('colina_data')

# ------------------------------------------------------------------------------

# load manager and loggers
mgr = PolymerManager(src_coll_path)
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER] # loggers from all modules which produce logging output

# filter for polymers with simulation and charges successfully generated
sim_poly = mgr.filtered_by(has_sims)
chg_methods = set(
    charge_method 
        for polymer in sim_poly.values()
            for charge_method in polymer.charges.keys() 
)

# create save folder and sets of structure converters and labels
save_dir.mkdir(exist_ok=True)
CHG_CVTRS = {
    'Charge heatmap 2D' : 'SMARTS',
    'Charge heatmap 3D' : 'CXSMARTS'
}

# iterate through simulated polymers and generate data and plots
for mol_name, polymer in sim_poly.items():
    # create folders to save to
    parent_dir = save_dir / mol_name
    parent_dir.mkdir(exist_ok=True)
    observ = defaultdict(defaultdict)

    fig_time, ax_time = None, None
    fig_acf, ax_acf = None, None

    # Plot relative charges
    for chg_pair in combinations(chg_methods, 2): # Plott all possible combinations of pairwise charge comparisons
        for label, cvtr in CHG_CVTRS.items():
            fig_chg, ax_chg = polymer.compare_charges(*sorted(chg_pair), converter=cvtr)
            fig_chg.savefig(parent_dir / f'{label} - {mol_name}', bbox_inches='tight')
            plt.close()
        
    # Plot observables
    for sim_dir, sim_paths_file in polymer.simulation_paths.items():
        sim_paths = SimulationPaths.from_file(sim_paths_file)
        if sim_paths.trajectory.suffix == '.dcd':
            # read simulation trajectory data
            sim_params = SimulationParameters.from_file(sim_paths.sim_params)
            chg_method = sim_params.charge_method
            time_data = pd.read_csv(sim_paths.time_data)
            spatial_data = pd.read_csv(sim_paths.spatial_data)

            # Plot RDFs
            radii, rdfs = trajectory.rdfs_to_plot_data(spatial_data)
            num_rdfs = len(rdfs.columns)
            nrows_rdf = num_rdfs // 3
            ncols_rdf = ceil(num_rdfs / nrows_rdf)
            fig_rdf, ax_rdf = plotutils.plot_df_props(radii, rdfs, nrows=nrows_rdf, ncols=ncols_rdf)
            fig_rdf.savefig(parent_dir / f'RDFs - {mol_name}.png', bbox_inches='tight')
            plt.close()

            # Plot time series' and autocorrelation functions
            times, time_series = trajectory.props_to_plot_data(time_data)
            if ax_time is None: # create new axes for first time...
                fig_time, ax_time = plotutils.plot_df_props(times, time_series)
            else: # otherwise, overlay plots
                for ax, (prop_name, data) in zip(ax_time, time_series.items()):
                    ax.plot(times, data)
                    ax.legend(chg_methods)

            time_series_acfs = time_series.apply(statistics.autocorrelate)
            time_series_acfs = time_series_acfs.rename(mapper=lambda name : re.sub('\(.*?\)', 'ACF', name), axis='columns')
            if ax_acf is None: # create new axes for first time...
                fig_acf, ax_acf = plotutils.plot_df_props(times, time_series_acfs)
            else: # otherwise, overlay plots
                for ax, (prop_name, data) in zip(ax_acf, time_series_acfs.items()):
                    ax.plot(times, data)
                    ax.legend(chg_methods)

            # Determine equilibration times and average observables
            nrows, ncols = 1, len(time_series.columns)
            fig_equil, ax_equil = plotutils.presize_subplots(nrows, ncols)
            for ax, (prop_name, series) in zip(ax_equil, time_series.items()):
                series = np.array(series)
                equil_idx = statistics.equil_loc(series)
                equil_time = times.to_numpy().flatten()[equil_idx]
                obs = round(np.mean(series[equil_idx:]), 3)
                observ[f'OpenFF 2.0.0 - {chg_method}'][prop_name] = obs

                ax.plot(times, series)
                ax.axvline(equil_time, color='r')
                ax.set_xlabel(times.columns[0])
                ax.set_ylabel(prop_name)
                ax.legend([f'<O> = {obs} (after t={equil_time})'], handlelength=0, loc='best') # position annotation in empty region
            fig_equil.savefig(parent_dir / f'Equil times {chg_method} - {mol_name}.png', bbox_inches='tight')
            plt.close()
            plt.close()

    # Plot overlaid property and ACF data
    fig_time.legend()
    fig_time.savefig(parent_dir / f'Polymer Properties - {mol_name}.png', bbox_inches='tight')
    plt.close()
    
    fig_acf.legend()
    fig_acf.savefig(parent_dir / f'Property ACFs - {mol_name}.png', bbox_inches='tight')
    plt.close()

    # Save CSV tables of observables, with experimental data appended
    compare_path = compare_dir / f'{mol_name}_exp.csv'
    observ_df = pd.DataFrame(observ)
    compare_df = pd.read_csv(compare_path, index_col=0) # gets rid of unnecessary index col, enables smooth concatenation

    observ_df_aug = pd.concat([observ_df, compare_df], axis=1)
    observ_df_aug.to_csv(parent_dir/ f'Observables - {mol_name}.csv')