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

def plot_normalised_densities(results):
    """
    Plot the normalised densities of elastin, collagen, immune cells, latent and active TGF Beta.
    """
    time = results["time"]

    sv_collagen_me = results["collagen_me"]
    sv_elastin_me = results["elastin_me"]
    sv_elastases = results["elastases"]
    sv_collagenases = results["collagenases"]
    sv_immune_cells = results["immune_cells"]
    sv_fibroblast = results["fibroblast"]
    sv_muscle_cells = results["muscle_cells"]
    sv_collagen_ad = results["collagen_ad"]
    sv_procollagen = results["procollagen"]
    sv_collagenase = results["collagenase"]
    sv_zymogen = results["zymogen"]
    sv_timp = results["timp"]
    sv_latent_tgf = results["latent_tgf_beta"]
    sv_active_tgf = results["active_tgf_beta"]

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.plot(time, sv_collagen_me, color='brown', linestyle=':', label='Medial Collagen', linewidth=1)
    ax.plot(time, sv_elastin_me, color='brown', linestyle=':', label='Medial Elastin', linewidth=1)
    ax.plot(time, sv_elastases, color='red', marker='o', markevery=300, label='Elastases', linewidth=1)
    ax.plot(time, sv_collagenases, color='red', marker='o', markevery=300, label='Collagenases', linewidth=1)
    ax.plot(time, sv_fibroblast, color='brown', label='Fibroblasts', linewidth=3)
    ax.plot(time, sv_muscle_cells, color='red', marker='v', markevery=300, label='Muscle Cells', linewidth=1)
    ax.plot(time, sv_collagen_ad, color='brown', label='Adventitial Collagen', linewidth=1)
    ax.plot(time, sv_procollagen, color='brown', label='Procollagen', linewidth=3)
    ax.plot(time, sv_collagenase, color='magenta', marker='^', markevery=300, linestyle='-', label='Collagenase', linewidth=1)
    ax.plot(time, sv_zymogen, color='magenta', marker='^', markevery=300, linestyle='--', label='Zymogen', linewidth=1)
    ax.plot(time, sv_timp, color='gold', linestyle='-', marker='v', markevery=300, label='TIMP', linewidth=1)
    ax.plot(time, sv_immune_cells, color='red', marker='o', markevery=300, label='Immune Cells', linewidth=3)
    ax.plot(time, sv_latent_tgf, color='green', linestyle='--', marker='s', markevery=300, label='Latent TGF-Beta', linewidth=1)
    ax.plot(time, sv_active_tgf, color='green', linestyle='-', marker='s', markevery=300, label='Active TGF-Beta', linewidth=1)

    ax.set_title('Normalised Densities Over Time', fontsize=16, weight='bold')
    ax.set_xlabel('Time (years)', fontsize=14)
    ax.set_ylabel('Normalised Density', fontsize=14)
    ax.axvline(40, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(50, color='black', linestyle=':', linewidth=1.5)
    ax.legend(loc='upper right', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim(38, 75)
    fig.tight_layout()
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
    

def plot_tgf_beta_treatment():
    """
    Placeholder for TGF-Beta treatment plot.
    """
    # This function is a placeholder for future implementation.
    # It should plot the effects of TGF-Beta treatment on arterial parameters.
    pass
