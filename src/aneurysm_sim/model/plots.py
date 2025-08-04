import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


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

def plot_normalised_densities(results, ax=None, title=None, legend=False, xlabel=False, ylabel=False):
    """
    Plot the normalised densities of elastin, collagen, immune cells, latent and active TGF Beta
    into the supplied Axes (or create one if none provided).
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
    Arrange 3 genotype plots in a 2x2 grid with the top-right
    subplot reserved for the legend box.
    """
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])

    ax1 = fig.add_subplot(gs[0, 0])  # Top-left: TT
    ax_legend = fig.add_subplot(gs[0, 1])  # Top-right: Legend placeholder
    ax2 = fig.add_subplot(gs[1, 0])  # Bottom-left: TC
    ax3 = fig.add_subplot(gs[1, 1])  # Bottom-right: CC

    # Plotting without legends on subplots
    plot_normalised_densities(results_tt, ax1, title="Genotype: TT", ylabel=True, xlabel=False, legend=False)
    plot_normalised_densities(results_tc, ax2, title="Genotype: TC", ylabel=True, xlabel=True, legend=False)
    plot_normalised_densities(results_cc, ax3, title="Genotype: CC", ylabel=False, xlabel=True, legend=False)

    # Remove axes for legend box
    ax_legend.axis('off')

    # Get legend handles/labels from one of the plots
    handles, labels = ax1.get_legend_handles_labels()

    # Add legend to the empty subplot area
    ax_legend.legend(handles, labels, loc='lower left', fontsize=10)

    # Shared axis labels for the figure
    # fig.text(0.5, 0.04, 'Time (years)', ha='center', fontsize=16)
    # fig.text(0.07, 0.5, 'Normalized Density', va='center', rotation='vertical', fontsize=16)

    plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.12, right=0.95, bottom=0.1, top=0.95)
    plt.show()

    return fig

    

def plot_systolic_stretch_over_time(results_tt, results_tc, results_cc):
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.plot(results_tt["time"], results_tt["lambda_sys"], label="TT", color='red')
    ax.plot(results_tc["time"], results_tc["lambda_sys"], label="TC", color='brown')
    ax.plot(results_cc["time"], results_cc["lambda_sys"], label="CC", color='gold')
    
    # ax.set_title("Systolic Stretch Over Time by Genotype (with treatment)", fontsize=16, weight='bold')
    ax.set_xlabel("Time (years)", fontsize=16)
    ax.set_ylabel(r"Systolic Stretch $\lambda_{sys}$", fontsize=16)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim(40, 75)
    ax.set_ylim(1.3, 1.5)
    ax.legend(fontsize=12)
    ax.axvline(45, color='black', linestyle='--', linewidth=1.5)
    fig.tight_layout()
    plt.show()
    
    return fig

def plot_max_collagen_stretch(results):
    time = results["time"]
    sv_lambda_c_max = results["lambd_c_max_history"]
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(time, sv_lambda_c_max, color='blue', label='Max Collagen Stretch', linewidth=2)
    # ax.set_title('Max Collagen Stretch Over Time', fontsize=16, weight='bold')
    ax.set_xlabel('Time (years)', fontsize=16)
    ax.set_ylabel('Max Collagen Stretch', fontsize=16)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim(40, 75)
    # ax.set_ylim(1.0, 1.2)
    ax.legend(fontsize=12)
    fig.tight_layout()
    plt.show()
    
    return fig

def plot_diameter_treatment_times(results_tt_list, results_tc_list, results_cc_list, t_treat_list):
    """
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    time = results_tt_list[0]["time"]
    for i, t_treat in enumerate(t_treat_list): 
        axes[0].plot(time, results_tt_list[i]["diameter"]/1e3, label=f"{t_treat} Years")
        axes[1].plot(time, results_tc_list[i]["diameter"]/1e3, label=f"{t_treat} Years")
        axes[2].plot(time, results_cc_list[i]["diameter"]/1e3, label=f"{t_treat} Years")

    axes[0].set_title('Genotype: TT')
    axes[1].set_title('Genotype: TC')
    axes[2].set_title('Genotype: CC')

    for ax in axes: 
        ax.set_xlabel('Time (years)', fontsize=16)
        ax.set_ylabel('Diameter (mm)', fontsize=16)
        ax.set_xlim(40, 75)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(title='TGF-β Treatment Time', fontsize=10)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title='t₀ (years)', loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.05), fontsize=11)
    plt.subplots_adjust(bottom=0.2, top=0.88)
    # fig.suptitle('Diameter Over Time with TGF-β Treatment', fontsize=16, weight='bold')
    plt.show()

    return fig

def plot_tgf_vs_fibroblasts(results_tt, results_tc, results_cc):
    fig, ax = plt.subplots(figsize=(7,7))

    # Package results for iteration
    all_results = {
        "TT": results_tt,
        "TC": results_tc,
        "CC": results_cc
    }

    for genotype, data in all_results.items():
        active_tgf = data['active_tgf_beta']
        fibroblasts = data['fibroblast']

        ax.plot(active_tgf, fibroblasts, label=genotype, lw=2)
        # mark start (o) and end (x) points
        ax.scatter(active_tgf[0], fibroblasts[0], marker="o", color=ax.lines[-1].get_color(), s=60)
        ax.scatter(active_tgf[-1], fibroblasts[-1], marker="x", color=ax.lines[-1].get_color(), s=80)

    ax.set_xlabel("Active TGF-β", fontsize=14)
    ax.set_ylabel("Fibroblasts", fontsize=14)
    # ax.set_title("Active TGF-β vs Fibroblasts (Phase Plane)", fontsize=16)
    ax.legend(title="Genotype")
    ax.grid(True)

    plt.show()

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
    treatment vs. no-treatment, for each genotype (TT, TC, CC), with a truly shared y-axis.
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

    fig, axes = plt.subplots(1, 3, figsize=(22, 7), sharey=True)

    for ax, aucs, title in zip(
        axes, [aucs_tt, aucs_tc, aucs_cc], ["Genotype: TT", "Genotype: TC", "Genotype: CC"]
    ):
        x = np.arange(len(components))
        width = 0.35

        ax.bar(x - width/2, aucs[:, 0], width, label=treatment_labels[0], color="skyblue")
        ax.bar(x + width/2, aucs[:, 1], width, label=treatment_labels[1], color="salmon")

        ax.set_title(title, fontsize=14, weight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([label for _, label in components], rotation=45, ha="right", fontsize=11)
        ax.grid(axis="y", linestyle="--", alpha=0.6)

    # Hide redundant y-axes (middle and right panels)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
        ax.spines['left'].set_visible(False)

    # Shared y-label for the entire figure
    fig.text(0.02, 0.5, "AUC", va="center", rotation="vertical", fontsize=14)

    axes[0].legend(loc="upper right", fontsize=11)
    plt.tight_layout(rect=[0.05, 0, 1, 1])  # leave space for shared ylabel
    plt.show()

    return fig




