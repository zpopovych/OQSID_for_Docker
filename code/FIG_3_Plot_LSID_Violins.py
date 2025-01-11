import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

γ_list = [ "0.079477",  "0.25133", "0.79477", "2.5133", "7.9477", "25.133", "79.477", "251.33"]

def read_fidelities(file):
    fidelities = []

    for γ_i in γ_list:

        with h5py.File(file, 'r') as f:
                
            fidelity_list = [f[γ_i][f"D{i}"]["Fidelity"][...] for i in range(1, 11)  if f[γ_i][f"D{i}"]["Fidelity"].size > 0]
            concatenated_list = [item for lst in fidelity_list for item in lst] 
            fidelities.append(concatenated_list) 
    
    return(fidelities)

def read_rank5_fidelities(file):
    fidelities_dmd = []
    fidelities_era = []

    for γ_i in γ_list:

        with h5py.File(file, 'r') as f:
                
            fidelity_dmd_list = [f[γ_i][f"D{i}"]["F_dmd_sb"][...] for i in range(1, 11)  if f[γ_i][f"D{i}"]["F_dmd_sb"].size > 0]
            concatenated_dmd_list = [item for lst in fidelity_dmd_list for item in lst] 
            fidelities_dmd.append(concatenated_dmd_list) 

            fidelity_era_list = [f[γ_i][f"D{i}"]["F_era_sb"][...] for i in range(1, 11)  if f[γ_i][f"D{i}"]["F_era_sb"].size > 0]
            concatenated_era_list = [item for lst in fidelity_era_list for item in lst] 
            fidelities_era.append(concatenated_era_list) 
    
    return(fidelities_dmd, fidelities_era)

dmd4_fids = read_fidelities("../results/DMD-Bloch4D_SB_trn4_tst10.h5")

fidelities_dmd5, fidelities_era5 = read_rank5_fidelities("../results/DMDvsERA_rank5_SB_trn4_tst10.h5")

with h5py.File("../results/MaxPosEigs_DMD_Bloch4.h5", "r") as file:
    PosEvals = file["PosEigs"][:]

import matplotlib.patches as mpatches

labels = []

def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))

Q = [316,  100,  31.6,   10,   3.16,    1,    0.316,    0.1]

labels = []

# Set the font sizes globally
plt.rcParams.update({'font.size': 14,  # Increase overall font size
                    'axes.labelsize': 16,  # Axis label size
                    'xtick.labelsize': 14,  # X-axis tick label size
                    'ytick.labelsize': 14,  # Y-axis tick label size
                    'legend.fontsize': 14,  # Legend font size
                    'axes.titlesize': 18})  # Title font size (if you have any)

# Create the first plot
fig, ax1 = plt.subplots(figsize=(8, 6))

add_label(ax1.violinplot([1-np.abs(_)  for _ in dmd4_fids], showextrema=True, widths=0.8), "(a) - DMD rank 4 \n (augmented Bloch)")

add_label(ax1.violinplot([1-np.abs(_)  for _ in fidelities_dmd5], showextrema=True, widths=0.8), "(b) - DMD rank 5 \n (SVD based)")

#showmedians=True, 

ax2 = ax1.twinx()
ax2.plot(np.arange(1,9), PosEvals, marker="+", color="red", label = "(c) - max{Re[eig($\mathcal{A}$)]} \n (DMD rank 4)")
#ax2.set_yscale('log')
ax2.set_ylim(0, 1e-10/2)
ax2.set_ylabel(r'Largest positive eigenvalue of $\mathcal{A}$ ', color='red')

ax1.set_yscale('log')
ax1.set_ylim(1e-6, 2.5)

#ax1.tight_layout()

ax1.set_xticks(range(1, len(γ_list) + 1), Q)
ax1.set_xlabel("$Q= \\nu / \\gamma$, quality factor")

ax1.set_ylabel("$|1 - F|$, infidelity")

ax1.legend(*zip(*labels), loc=2)

# Changing the color of the second y-axis
ax2.tick_params(axis='y', colors='red')
ax2.spines['right'].set_color('red')
ax2.legend(loc=1)

#ax1.set_title("Performance of Markovian models for spin boson system")

fig.savefig("../results/FIG_3_DMD_Bloch4-5.pdf")

print("Ploting Linear SID violins done.")