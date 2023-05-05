'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
import re
from pathlib import Path
from collections import defaultdict

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

# ------------------------------------------------------------------------------

mgr = PolymerManager(src_coll_path)
main_logger = logging.getLogger(__name__)
loggers = [main_logger, *LOGGERS_MASTER] # loggers from all modules which produce logging output

save_dir.mkdir(exist_ok=True)
heatmap_dir = save_dir / 'heatmaps'
heatmap_dir.mkdir(exist_ok=True)

# filter for polymers with simulation and charges successfully generated
sim_poly = mgr.filtered_by(has_sims)
chg_methods = set(
    charge_method 
        for polymer in sim_poly.values()
            for charge_method in polymer.charges.keys() 
)

# create workspaces for each charge method
observ_sets, chg_dirs = {}, {}
for chg_method in chg_methods:
    observ_sets[chg_method] = defaultdict(defaultdict) # to accumulate averaged observables
    
    chg_dir = save_dir / chg_method # to save plots and data to
    chg_dir.mkdir(exist_ok=True)
    chg_dirs[chg_method] = chg_dir

# iterate through simulated polymers and generate data and plots
for mol_name, polymer in sim_poly.items():
    for sim_dir, sim_paths_file in polymer.simulation_paths.items():
        sim_paths = SimulationPaths.from_file(sim_paths_file)
        if sim_paths.trajectory.suffix == '.dcd':
            # read simulation trajectory data
            sim_params = SimulationParameters.from_file(sim_paths.sim_params)
            time_data = pd.read_csv(sim_paths.time_data)
            spatial_data = pd.read_csv(sim_paths.spatial_data)

            # create folders to save to
            parent_dir = chg_dirs[sim_params.charge_method] / mol_name
            parent_dir.mkdir(exist_ok=True)

            # Plot relative charges
            fig_chg_2d, ax_chg_2d = polymer.compare_charges('ABE10_averaged', 'Espaloma_AM1BCC', converter='SMARTS')
            fig_chg_2d.savefig(heatmap_dir / f'Charge heatmap 2D - {mol_name}', bbox_inches='tight')
            plt.close()
            
            fig_chg_3d, ax_chg_3d = polymer.compare_charges('ABE10_averaged', 'Espaloma_AM1BCC', converter='CXSMARTS')
            fig_chg_3d.savefig(heatmap_dir / f'Charge heatmap 3D - {mol_name}', bbox_inches='tight')
            plt.close()

            # Plot RDFs
            radii, rdfs = trajectory.rdfs_to_plot_data(spatial_data)
            num_rdfs = len(rdfs.columns)
            nrows_rdf = num_rdfs // 3
            ncols_rdf = ceil(num_rdfs / nrows_rdf)
            fig_rdf, ax_rdf = plotutils.plot_df_props(radii, rdfs, nrows=nrows_rdf, ncols=ncols_rdf)
            fig_rdf.savefig(parent_dir / f'RDFs - {mol_name}.png', bbox_inches='tight')
            plt.close()

            # For storing a table of averaged observables
            observ = defaultdict(defaultdict)

            # Plot time series' and autocorrelation functions
            times, time_series = trajectory.props_to_plot_data(time_data)
            fig_time, ax_time = plotutils.plot_df_props(times, time_series)
            fig_time.savefig(parent_dir / f'Polymer Properties - {mol_name}.png', bbox_inches='tight')
            plt.close()

            time_series_acfs = time_series.apply(statistics.autocorrelate)
            time_series_acfs = time_series_acfs.rename(mapper=lambda name : re.sub('\(.*?\)', 'ACF', name), axis='columns')
            fig_acf, ax_acf = plotutils.plot_df_props(times, time_series_acfs)
            fig_acf.savefig(parent_dir / f'Property ACFs - {mol_name}.png', bbox_inches='tight')
            plt.close()

            # Determine equilibration times and average observables
            nrows, ncols = 1, len(time_series.columns)
            fig_equil, ax_equil = plotutils.presize_subplots(nrows, ncols)
            for ax, (prop_name, series) in zip(ax_equil, time_series.items()):
                series = np.array(series)
                equil_idx = statistics.equil_loc(series)
                equil_time = times.to_numpy().flatten()[equil_idx]
                obs = round(np.mean(series[equil_idx:]), 3)
                observ_sets[sim_params.charge_method][polymer.base_mol_name][prop_name] = obs

                ax.plot(times, series)
                ax.axvline(equil_time, color='r')
                ax.set_xlabel(times.columns[0])
                ax.set_ylabel(prop_name)
                ax.legend([f'<O> = {obs} (after t={equil_time})'], handlelength=0, loc='best') # position annotation in empty region
            fig_equil.savefig(parent_dir / f'Equilibration times - {mol_name}.png', bbox_inches='tight')
            plt.close()
            plt.close()

# save tables of observables
for chg_method, observ in observ_sets.items():
    chg_dir = chg_dirs[chg_method]
    observ_df = pd.DataFrame(observ)
    observ_df.to_csv(chg_dir / f'Observables - {chg_method}.csv')