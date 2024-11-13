import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-poster')

# Connect to SQLite database
db_path = './data/xy_model.db'
conn = sqlite3.connect(db_path)

# Load data into a pandas DataFrame
query = "SELECT * FROM SimulationResults"
df = pd.read_sql_query(query, conn)

# Close the database connection
conn.close()

# Calculate averages for energy, magnetization, and susceptibility per temperature and lattice size
df_avg = df.groupby(['T', 'L']).mean().reset_index()

# Extract unique temperature and lattice sizes
T_range = sorted(df_avg['T'].unique())
lattice_sizes = sorted(df_avg['L'].unique())

# Initialize dictionaries to hold data for different lattice sizes
energy_data = {}
magnetization_data = {}
susceptibility_data = {}
vortex_data = {}
# Populate dictionaries with averaged data
for i, L in enumerate(lattice_sizes):
    energy_data[L] = df_avg[df_avg['L'] == L]['energy'].values
    magnetization_data[L] = df_avg[df_avg['L'] == L]['magnetization'].values
    susceptibility_data[L] = df_avg[df_avg['L'] == L]['susceptibility'].values
    vortex_data[L] = df_avg[df_avg['L'] == L]['vortex_count'].values

label = ['o','^','s','P','D','2','*','p']

# Plot energy vs. temperature
plt.figure(figsize=(12, 8))
for i, L in enumerate(lattice_sizes):
    plt.plot(T_range, energy_data[L], label=f'L={L}', marker=f'{label[i]}', markerfacecolor='none', markeredgewidth=1.0, linestyle='none')
plt.xlabel("$T$", fontsize=16)
plt.ylabel("$E$", fontsize=16)
plt.legend(loc='upper left')
plt.grid(True, linewidth=1, color='#c7c7c7')
plt.savefig("plots/energy_v_t.png", dpi=500, bbox_inches='tight')
plt.show()

# Plot magnetization vs. temperature
plt.figure(figsize=(12, 8))
for i, L in enumerate(lattice_sizes):
    plt.plot(T_range, magnetization_data[L], label=f'L={L}', marker=f'{label[i]}', markerfacecolor='none', markeredgewidth=1.0, linestyle='none')
plt.xlabel("$T$", fontsize=16)
plt.ylabel("$M$", fontsize=16)
plt.grid(True, linewidth=1, color='#c7c7c7')
plt.legend()
plt.savefig("plots/mag_v_t.png", dpi=500, bbox_inches='tight')
plt.show()

# Plot susceptibility vs. temperature (log scale)
plt.figure(figsize=(12, 8))
for i, L in enumerate(lattice_sizes):
    plt.plot(T_range, susceptibility_data[L], label=f'L={L}', marker=f'{label[i]}', markerfacecolor='none', markeredgewidth=1.0, linestyle='none')
plt.xlabel("$T$", fontsize=16)
plt.ylabel("$S$", fontsize=16)
plt.grid(True, linewidth=1, color='#c7c7c7')
plt.legend()
plt.yscale('log')
plt.savefig("plots/sus_v_t.png", dpi=500, bbox_inches='tight')
plt.show()

# Plot vortex count vs. temperature
plt.figure(figsize=(12, 8))
for i, L in enumerate(lattice_sizes):
    plt.plot(T_range, vortex_data[L], label=f'L={L}', marker=f'{label[i]}', markerfacecolor='none', markeredgewidth=1.0, linestyle='none')
plt.xlabel("$T$", fontsize=16)
plt.ylabel("Vortex Count", fontsize=16)
plt.grid(True, linewidth=1, color='#c7c7c7')
plt.legend()
plt.savefig("plots/vortex_v_t.png", dpi=500, bbox_inches='tight')
plt.show()