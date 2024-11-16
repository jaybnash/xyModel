import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scienceplots

# Python script to plot data from xy model simulations
# Accesses database to retrieve simulation results
# Uses matplotlib to create plots of energy, magnetization, susceptibility, and vortex count vs. temperature
# Since the raw data is not submitted, a sample plot is provided in the submission files
# Jay Nash 11/15/2024

# Style the plots so they look nice
plt.style.use(['science','no-latex'])

# Connect to SQLite database
db_path = './data/xy_model.db'
conn = sqlite3.connect(db_path)

# Load data into a pandas DataFrame
query = "SELECT * FROM SimulationResults"
df = pd.read_sql_query(query, conn)

# Close the database connection
conn.close()

# Averages energy, magnetization, susceptibility, and vortex count per temperature and lattice size
df_avg = df.groupby(['T', 'L']).mean().reset_index()

# Get unique temperature and lattice sizes
T_range = sorted(df_avg['T'].unique())
lattice_sizes = sorted(df_avg['L'].unique())

# Init dictionaries to hold data for different lattice sizes
energy_data = {}
magnetization_data = {}
susceptibility_data = {}
vortex_data = {}

# Make dictionaries with averaged data
for i, L in enumerate(lattice_sizes):
    energy_data[L] = df_avg[df_avg['L'] == L]['energy'].values
    magnetization_data[L] = df_avg[df_avg['L'] == L]['magnetization'].values
    susceptibility_data[L] = df_avg[df_avg['L'] == L]['susceptibility'].values
    vortex_data[L] = df_avg[df_avg['L'] == L]['vortex_count'].values

label = ['o','^','s','P','D','2','*','p']

# Create a 2x2 subplot grid
fig, axs = plt.subplots(2, 2, figsize=(12, 12))

# Plot energy vs. temperature
for i, L in enumerate(lattice_sizes):
    axs[0, 0].plot(T_range, energy_data[L], label=f'L={L}', marker=f'{label[i]}', markerfacecolor='none', markeredgewidth=1.0, linestyle='none')
axs[0, 0].set_xlabel("$T$", fontsize=16)
axs[0, 0].set_ylabel("$E$", fontsize=16)
axs[0, 0].legend(title='Lattice Size', title_fontsize=12)
axs[0, 0].grid(True, linewidth=1, color='#c7c7c7')
axs[0, 0].set_title("Energy vs. Temperature")

# Plot magnetization vs. temperature
for i, L in enumerate(lattice_sizes):
    axs[0, 1].plot(T_range, magnetization_data[L], label=f'L={L}', marker=f'{label[i]}', markerfacecolor='none', markeredgewidth=1.0, linestyle='none')
axs[0, 1].set_xlabel("$T$", fontsize=16)
axs[0, 1].set_ylabel("$M$", fontsize=16)
axs[0, 1].grid(True, linewidth=1, color='#c7c7c7')
axs[0, 1].legend(title='Lattice Size', title_fontsize=12)
axs[0, 1].set_title("Magnetization vs. Temperature")

# Plot susceptibility vs. temperature (log scale)
for i, L in enumerate(lattice_sizes):
    axs[1, 0].plot(T_range, susceptibility_data[L], label=f'L={L}', marker=f'{label[i]}', markerfacecolor='none', markeredgewidth=1.0, linestyle='none')
axs[1, 0].set_xlabel("$T$", fontsize=16)
axs[1, 0].set_ylabel("$S$", fontsize=16)
axs[1, 0].grid(True, linewidth=1, color='#c7c7c7')
axs[1, 0].legend(title='Lattice Size', title_fontsize=12)
axs[1, 0].set_yscale('log')
axs[1, 0].set_title("Susceptibility vs. Temperature")

# Plot vortex count vs. temperature
for i, L in enumerate(lattice_sizes):
    axs[1, 1].plot(T_range, vortex_data[L], label=f'L={L}', marker=f'{label[i]}', markerfacecolor='none', markeredgewidth=1.0, linestyle='none')
axs[1, 1].set_xlabel("$T$", fontsize=16)
axs[1, 1].set_ylabel("$V$", fontsize=16)
axs[1, 1].grid(True, linewidth=1, color='#c7c7c7')
axs[1, 1].legend(title='Lattice Size', title_fontsize=12)
axs[1, 1].set_title("Vortex Count vs. Temperature")

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig("plots/combined_plots.png", dpi=500, bbox_inches='tight')
plt.show()
