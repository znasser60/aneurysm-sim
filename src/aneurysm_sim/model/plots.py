import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def plot_pressure_vs_stretch(results, n_zoom=120): 
    """
    Plot the pressure vs stretch of the artery with a zoomed-in view.
    """

    sv_stretch_var = results["sv_stretch_var"]
    sv_pressure_var = results["sv_pressure_var"]
    sv_pressure_var_elastin = results["sv_pressure_var_elastin"]
    sv_pressure_var_collagen = results["sv_pressure_var_collagen"]
    sv_pressure_var_muscle_p = results["sv_pressure_var_muscle_p"]
    sv_pressure_var_muscle_a = results["sv_pressure_var_muscle_a"]

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
    Plot the pressure vs diameter of the artery with a zoomed-in view.
    """

    sv_diam_var = results["sv_diam_var"]
    sv_pressure_var = results["sv_pressure_var"]
    sv_pressure_var_elastin = results["sv_pressure_var_elastin"]
    sv_pressure_var_collagen = results["sv_pressure_var_collagen"]
    sv_pressure_var_muscle_p = results["sv_pressure_var_muscle_p"]
    sv_pressure_var_muscle_a = results["sv_pressure_var_muscle_a"]


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

def plot_stretch_vs_stress(results, n_zoom=120):
    """
    Plot the stretch vs stress with a zoomed-in view.
    """

    sv_stretch_var = results["sv_stretch_var"]
    sv_stress_var_elastin = results["sv_stress_var_elastin"]
    sv_stress_var_collagen_me = results["sv_stress_var_collagen_me"]
    sv_stress_var_collagen_ad = results["sv_stress_var_collagen_ad"]
    sv_stress_var_collagen = results["sv_stress_var_collagen"]
    sv_stress_var_total = results["sv_stress_var_total"]

    plt.figure(figsize=(12, 8))
    plt.plot(sv_stretch_var[32:n_zoom], sv_stress_var_elastin[32:n_zoom], '--', linewidth=2, label='Elastin')
    plt.plot(sv_stretch_var[64:n_zoom], sv_stress_var_collagen[64:n_zoom], '-.', linewidth=2, label='Collagen Total')
    plt.plot(sv_stretch_var[44:n_zoom], sv_stress_var_collagen_me[44:n_zoom], '+', linewidth=2, markersize=8, label='Collagen ME')
    plt.plot(sv_stretch_var[:n_zoom], sv_stress_var_collagen_ad[:n_zoom], '.', linewidth=2, markersize=8, label='Collagen AD')
    plt.plot(sv_stretch_var[:n_zoom], sv_stress_var_total[:n_zoom], '-', linewidth=3, label='Total Stress')
    plt.xlabel(r'Stretch $\lambda$', fontsize=14)
    plt.ylabel(r'Stress $\sigma$', fontsize=14)
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.show()

def plot_normalised_densities(results, ax=None, title=None, legend=False, xlabel=False, ylabel=False):
    """
    Plot the normalised densities of elastin, collagen, immune cells, latent and active TGF Beta.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        xlabel, ylabel = True, True


    time = results["time"]

    ax.plot(time, results["collagen_me"],   color='brown',   linestyle=':', label='Medial Collagen',    linewidth=1)
    ax.plot(time, results["elastin_me"],    color='brown',   linestyle=':', label='Medial Elastin',     linewidth=1)
    ax.plot(time, results["elastases"],     color='red',     marker='o',  markevery=300, label='Elastases',       linewidth=1)
    ax.plot(time, results["collagenases"],  color='red',     marker='o',  markevery=300, label='Collagenases',    linewidth=1)
    ax.plot(time, results["fibroblast"],    color='brown',   linewidth=3, label='Fibroblasts')
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
    
    if xlabel:
        ax.set_xlabel('Time (years)', fontsize=14)
    if ylabel:
        ax.set_ylabel('Normalized Density', fontsize=14)

    ax.set_xlim(38, 75)
    ax.axvline(40, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(50, color='black', linestyle=':',  linewidth=1.5)
    ax.grid(True, linestyle='--', alpha=0.4)
    if legend:
        ax.legend(loc='upper right', fontsize=10, ncol=2)
    return ax

def plot_normalised_densities_by_genotype(results_tt, results_tc, results_cc):
    """
    Plot normalised densities for each genotype (TT, TC, CC) in a 2x2 grid.
    """
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])

    ax1 = fig.add_subplot(gs[0, 0])  # Top-left: TT
    ax_legend = fig.add_subplot(gs[0, 1])  # Top-right: Legend placeholder
    ax2 = fig.add_subplot(gs[1, 0])  # Bottom-left: TC
    ax3 = fig.add_subplot(gs[1, 1])  # Bottom-right: CC

    plot_normalised_densities(results_tt, ax1, title="Genotype: TT", ylabel=True, xlabel=False, legend=False)
    plot_normalised_densities(results_tc, ax2, title="Genotype: TC", ylabel=True, xlabel=True, legend=False)
    plot_normalised_densities(results_cc, ax3, title="Genotype: CC", ylabel=False, xlabel=True, legend=False)
    ax_legend.axis('off')
    handles, labels = ax1.get_legend_handles_labels()
    ax_legend.legend(handles, labels, loc='lower left', fontsize=10)

    plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.12, right=0.95, bottom=0.1, top=0.95)
    plt.show()

    return fig

def plot_systolic_stretch_over_time(results_tt, results_tc, results_cc):
    """
    Plot the systolic stretch over time for each genotype.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.plot(results_tt["time"], results_tt["lambda_sys"], label="TT", color='red')
    ax.plot(results_tc["time"], results_tc["lambda_sys"], label="TC", color='brown')
    ax.plot(results_cc["time"], results_cc["lambda_sys"], label="CC", color='gold')
    
    # ax.set_title("Systolic Stretch Over Time by Genotype (with treatment)", fontsize=16, weight='bold')
    ax.set_xlabel("Time (years)", fontsize=20)
    ax.set_ylabel(r"Systolic Stretch $\lambda_{sys}$", fontsize=20)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_xlim(40, 75)
    # ax.set_ylim(1.3, 1.5)
    ax.legend(fontsize=12)
    ax.axvline(45, color='black', linestyle='--', linewidth=1.5)
    fig.tight_layout()
    plt.show()
    
    return fig

def compute_auc(time, y):
    """Compute area under the curve using trapezoidal rule."""
    return np.trapz(y, time)

def extract_aucs(results_dict, components, treatment_labels):
    """Return AUC values for each component and treatment condition."""
    aucs = []
    for var_key, _ in components:
        auc_vals = []
        for cond in treatment_labels:
            time = results_dict[cond]["time"]
            y = results_dict[cond][var_key]
            auc_vals.append(compute_auc(time, y))
        aucs.append(auc_vals)
    return np.array(aucs)

def plot_auc_bars_by_genotype(results_tt, results_tc, results_cc, treatment_labels=("No Treatment", "Treatment")):
    """
    Plot bar charts of AUC values for selected normalised densities, comparing
    treatment vs. no-treatment, for each genotype (TT, TC, CC).
    """

    components = [
        ("fibroblast", "Fibroblasts"),
        ("collagen_ad", "Adventitial Collagen"),
        ("latent_tgf_beta", "Latent TGF-β"),
        ("active_tgf_beta", "Active TGF-β"),
    ]

    aucs_tt = extract_aucs(results_tt, components, treatment_labels)
    aucs_tc = extract_aucs(results_tc, components, treatment_labels)
    aucs_cc = extract_aucs(results_cc, components, treatment_labels)
    print("AUCs (TT):", aucs_tt)
    print("AUCs (TC):", aucs_tc)
    print("AUCs (CC):", aucs_cc)

    fig, axes = plt.subplots(1, 3, figsize=(22, 7), sharey=True)

    for ax, aucs, title in zip(
        axes, [aucs_tt, aucs_tc, aucs_cc], ["Genotype: TT", "Genotype: TC", "Genotype: CC"]
    ):
        x = np.arange(len(components))
        width = 0.35
        bars1 = ax.bar(x - width/2, aucs[:, 0], width, label=treatment_labels[0], color="skyblue")
        bars2 = ax.bar(x + width/2, aucs[:, 1], width, label=treatment_labels[1], color="salmon")

        ax.set_title(title, fontsize=14, weight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([label for _, label in components], rotation=45, ha="right", fontsize=11)
        ax.grid(axis="y", linestyle="--", alpha=0.6)

    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
        ax.spines['left'].set_visible(False)

    fig.text(0.02, 0.5, "AUC", va="center", rotation="vertical", fontsize=14)

    axes[0].legend(loc="upper right", fontsize=11)
    plt.tight_layout(rect=[0.05, 0, 1, 1])  
    plt.show()

    return fig

def plot_att_dist(lambda_att_min, lambda_att_mode, lambda_att_max): 
    """
    Plot triangular distribution of attachment stretch.
    """
    x = np.linspace(lambda_att_min, lambda_att_max, 1000)
    y = np.where(
        (x >= lambda_att_min) & (x <= lambda_att_mode),
        2 * (x - lambda_att_min) / ((lambda_att_mode - lambda_att_min) * (lambda_att_max - lambda_att_min)),
        np.where(
            (x > lambda_att_mode) & (x <= lambda_att_max),
            2 * (lambda_att_max - x) / ((lambda_att_max - lambda_att_mode) * (lambda_att_max - lambda_att_min)),
            0
        )
    )

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='Attachment Stretch Distribution', color='blue')
    plt.xlabel(r'Attachment Stretch $\lambda_A^{AT}$', fontsize=14)
    plt.ylabel('PDF', fontsize=14)
    plt.grid(True)
    plt.show()

def plot_rec_dist(lambda_rec_min, lambda_rec_mode, lambda_rec_max): 
    """
    Plot triangular distribution of recruitment stretch.
    """
    x = np.linspace(lambda_rec_min, lambda_rec_max, 1000)
    y = np.where(
        (x >= lambda_rec_min) & (x <= lambda_rec_mode),
        2 * (x - lambda_rec_min) / ((lambda_rec_mode - lambda_rec_min) * (lambda_rec_max - lambda_rec_min)),
        np.where(
            (x > lambda_rec_mode) & (x <= lambda_rec_max),
            2 * (lambda_rec_max - x) / ((lambda_rec_max - lambda_rec_mode) * (lambda_rec_max - lambda_rec_min)),
            0
        )
    )

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='Recruitment Stretch Distribution', color='blue')
    plt.xlabel(r'Recruitment Stretch $\lambda_A^R$', fontsize=14)
    plt.ylabel('PDF', fontsize=14)
    plt.grid(True)
    plt.show()





