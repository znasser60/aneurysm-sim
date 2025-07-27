import numpy as np
import matplotlib.pyplot as plt


def plot_pressure_vs_stretch(results, n_zoom=120): 
    """
    Plot the pressure vs stretch with a zoomed-in view.
    """

    sv_stretch_var = results["sv_stretch_var"]
    sv_pressure_var = results["sv_pressure_var"]
    sv_pressure_var_elastin = results["sv_pressure_var_elastin"]
    sv_pressure_var_collagen = results["sv_pressure_var_collagen"]
    sv_pressure_var_collagen_me = results["sv_pressure_var_collagen_me"]
    sv_pressure_var_collagen_ad = results["sv_pressure_var_collagen_ad"]
    sv_pressure_var_muscle_p = results["sv_pressure_var_muscle_p"]
    sv_pressure_var_muscle_a = results["sv_pressure_var_muscle_a"]
    sv_pressure_var_coll = sv_pressure_var_collagen_me + sv_pressure_var_collagen_ad

    plt.figure(figsize=(12, 8))
    plt.plot(sv_stretch_var[:n_zoom], sv_pressure_var[:n_zoom]/1e3, linewidth=2, label='Total')
    plt.plot(sv_stretch_var[32:n_zoom], sv_pressure_var_elastin[32:n_zoom]/1e3, '--', linewidth=2, label='Elastin')
    plt.plot(sv_stretch_var[64:n_zoom], sv_pressure_var_collagen[64:n_zoom]/1e3, '-.', linewidth=2, label='Collagen')
    plt.plot(sv_stretch_var[44:n_zoom], sv_pressure_var_muscle_p[44:n_zoom]/1e3, '--', linewidth=2, label='Muscle Passive')
    plt.plot(sv_stretch_var[:n_zoom], sv_pressure_var_muscle_a[:n_zoom]/1e3, '--', linewidth=2, label='Muscle Active')
    plt.title('Pressure vs Stretch')
    plt.xlabel('Stretch')
    plt.ylabel('Pressure (kPa)')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.ylim([0, 60])
    plt.show()

def plot_pressure_vs_diameter(results, n_zoom=120): 
    """
    Plot the pressure vs diameter with a zoomed-in view."""

    sv_diam_var = results["sv_diam_var"]
    sv_pressure_var = results["sv_pressure_var"]
    sv_pressure_var_elastin = results["sv_pressure_var_elastin"]
    sv_pressure_var_collagen_me = results["sv_pressure_var_collagen_me"]
    sv_pressure_var_collagen_ad = results["sv_pressure_var_collagen_ad"]
    sv_pressure_var_collagen = results["sv_pressure_var_collagen"]
    sv_pressure_var_muscle_p = results["sv_pressure_var_muscle_p"]
    sv_pressure_var_muscle_a = results["sv_pressure_var_muscle_a"]

    sv_pressure_var_coll = sv_pressure_var_collagen_me + sv_pressure_var_collagen_ad

    plt.figure(figsize=(12, 8))
    plt.plot(sv_diam_var[:n_zoom], sv_pressure_var[:n_zoom]/1e3, '-', linewidth=6, label='Total')
    plt.plot(sv_diam_var[32:n_zoom], sv_pressure_var_elastin[32:n_zoom]/1e3, '--', linewidth=3, label='E')
    plt.plot(sv_diam_var[64:n_zoom], sv_pressure_var_collagen[64:n_zoom]/1e3, '-.', linewidth=3, label='C')
    plt.plot(sv_diam_var[44:n_zoom], sv_pressure_var_muscle_p[44:n_zoom]/1e3, '+', linewidth=2, markersize=8, label='VSMCp')
    plt.plot(sv_diam_var[:n_zoom], sv_pressure_var_muscle_a[:n_zoom]/1e3, '.', linewidth=2, markersize=8, label='VSMCa')

    plt.axvline(x=2.9, color='red', linestyle='--', linewidth=2)
    plt.axhline(y=16, color='red', linestyle='--', linewidth=2)
    plt.plot(2.9, 16, 'ro', linewidth=3, markersize=8)

    plt.title('Pressure vs Diameter')
    plt.xlabel('Diameter (mm)')
    plt.ylabel('Pressure (kPa)')
    plt.legend(loc='upper left')
    plt.ylim([0, 60])
    plt.grid(True)
    plt.show()

def plot_elastin_degradation(results_dict, n_zoom=120):
    """
    Plot pressure-stretch curves for different elastin degradation levels.
    """
    plt.figure(figsize=(10, 6))
    
    for label, results in results_dict.items():
        stretch = results["sv_stretch_var"][:n_zoom]
        pressure = results["sv_pressure_var"][:n_zoom] / 1e3  
        plt.plot(stretch, pressure, label=label, alpha=0.9)

    plt.title('Pressure vs. Stretch — Elastin Degradation', fontsize=16, weight='bold')
    plt.xlabel('Stretch', fontsize=14)
    plt.ylabel('Pressure (kPa)', fontsize=14)
    plt.ylim(0, 60)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend(title='Degradation Level', loc='upper left', fontsize=11)
    plt.tight_layout()
    plt.show()

def plot_normalised_densities(results, ax=None, title=None, legend=False):
    """
    Plot the normalised densities of elastin, collagen, immune cells, latent and active TGF Beta
    into the supplied Axes (or create one if none provided).
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    time = results["time"]

    ax.plot(time, results["collagen_me"],   color='brown',   linestyle=':', label='Medial Collagen',    linewidth=1)
    ax.plot(time, results["elastin_me"],    color='brown',   linestyle=':', label='Medial Elastin',     linewidth=1)
    ax.plot(time, results["elastases"],     color='red',     marker='o',  markevery=300, label='Elastases',       linewidth=1)
    ax.plot(time, results["collagenases"],  color='red',     marker='o',  markevery=300, label='Collagenases',    linewidth=1)
    ax.plot(time, results["fibroblast"],    color='brown',   linewidth=3, label='Fibroblasts')
    # ax.plot(time, results["muscle_cells"],  color='red',     marker='v',  markevery=300, label='Muscle Cells',    linewidth=1)
    ax.plot(time, results["collagen_ad"],   color='brown',   linewidth=1, label='Adventitial Collagen')
    ax.plot(time, results["procollagen"],   color='brown',   linewidth=3, label='Procollagen')
    ax.plot(time, results["collagenase"],   color='magenta', marker='^',  markevery=300, linestyle='-',  label='Collagenase',linewidth=1)
    ax.plot(time, results["zymogen"],       color='magenta', marker='^',  markevery=300, linestyle='--', label='Zymogen',     linewidth=1)
    ax.plot(time, results["timp"],          color='gold',    marker='v',  markevery=300, linestyle='-',  label='TIMP',        linewidth=1)
    ax.plot(time, results["immune_cells"],  color='red',     marker='o',  markevery=300, label='Immune Cells',   linewidth=3)
    ax.plot(time, results["latent_tgf_beta"], color='green', linestyle='--', marker='s', markevery=300, label='Latent TGF-β',linewidth=1)
    ax.plot(time, results["active_tgf_beta"], color='green', linestyle='-',  marker='s', markevery=300, label='Active TGF-β',linewidth=1)

    if title:
        ax.set_title(title, fontsize=14, weight='bold')

    ax.set_xlabel('Time (years)', fontsize=12)
    ax.set_ylabel('Normalized Density', fontsize=12)
    ax.set_xlim(38, 75)
    ax.axvline(40, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(50, color='black', linestyle=':',  linewidth=1.5)
    ax.grid(True, linestyle='--', alpha=0.4)
    if legend:
        ax.legend(loc='upper right', fontsize=10, ncol=2)
    return ax

def plot_normalised_densities_by_genotype(results_tt, results_tc, results_cc):
    """
    Compare TT, TC, and CC genotypes side by side in a single row.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    plot_normalised_densities(results_tt, axes[0], title="Genotype: TT")
    plot_normalised_densities(results_tc, axes[1], title="Genotype: TC")
    plot_normalised_densities(results_cc, axes[2], title="Genotype: CC")

    # axes[1].set_xlabel("Time (years)", fontsize=12)
    # axes[0].set_ylabel("Normalised Density", fontsize=12)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.01), fontsize=11)
    plt.subplots_adjust(bottom=0.25, top=0.9)
    fig.suptitle("Normalised Densities by TGF-β Genotype", fontsize=16, weight='bold')
    # fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    return fig
    

def plot_systolic_stretch_over_time(results_tt, results_tc, results_cc):
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.plot(results_tt["time"], results_tt["lambda_sys"], label="TT", color='red')
    ax.plot(results_tc["time"], results_tc["lambda_sys"], label="TC", color='brown')
    ax.plot(results_cc["time"], results_cc["lambda_sys"], label="CC", color='gold')
    
    ax.set_title("Systolic Stretch Over Time by Genotype", fontsize=16, weight='bold')
    ax.set_xlabel("Time (years)", fontsize=14)
    ax.set_ylabel("Systolic Stretch", fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim(40, 75)
    ax.set_ylim(1.3, 1.5)
    ax.legend(fontsize=12)
    fig.tight_layout()
    plt.show()
    
    return fig

def plot_max_collagen_stretch(results):
    time = results["time"]
    sv_lambda_c_max = results["lambd_c_max_history"]
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(time, sv_lambda_c_max, color='blue', label='Max Collagen Stretch', linewidth=2)
    ax.set_title('Max Collagen Stretch Over Time', fontsize=16, weight='bold')
    ax.set_xlabel('Time (years)', fontsize=14)
    ax.set_ylabel('Max Collagen Stretch', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim(40, 75)
    # ax.set_ylim(1.0, 1.2)
    ax.legend(fontsize=12)
    fig.tight_layout()
    plt.show()
    
    return fig
    
