import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

from matplotlib.colors import TwoSlopeNorm
from aneurysm_sim.model import functions
from aneurysm_sim.config.parameters import ArterialParameters

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
    plt.plot(sv_diam_var[32:n_zoom], sv_pressure_var_elastin[32:n_zoom]/1e3, '--', linewidth=3, label='Elastin')
    plt.plot(sv_diam_var[64:n_zoom], sv_pressure_var_collagen[64:n_zoom]/1e3, '-.', linewidth=3, label='Collagen')
    plt.plot(sv_diam_var[44:n_zoom], sv_pressure_var_muscle_p[44:n_zoom]/1e3, '+', linewidth=2, markersize=8, label='VSMCp')
    plt.plot(sv_diam_var[:n_zoom], sv_pressure_var_muscle_a[:n_zoom]/1e3, '.', linewidth=2, markersize=8, label='VSMCa')

    plt.axvline(x=2.9 , color='red', linestyle='--', linewidth=2)
    plt.axhline(y=16, color='red', linestyle='--', linewidth=2)
    plt.plot(2.9, 16, 'ro', linewidth=3, markersize=8)

    plt.title('Pressure vs Diameter', fontsize=30, weight='bold', pad=25)
    plt.xlabel('Diameter (mm)', fontsize=20, weight='bold', labelpad=15)
    plt.ylabel('Pressure (kPa)', fontsize=20, weight='bold', labelpad=15)
    plt.legend(loc='upper left', fontsize=16, frameon=True, shadow=True, edgecolor='black', title='Stress Components', title_fontsize=18)
    plt.ylim([0, 60])
    plt.tick_params(axis='both', which='major', labelsize=18)
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
    sv_stress_var_muscle = results["sv_stress_var_muscle_t"]
    sv_stress_var_total = results["sv_stress_var_total"]

    plt.figure(figsize=(12, 8))
    plt.plot(sv_stretch_var[32:n_zoom], sv_stress_var_elastin[32:n_zoom], '--', linewidth=2, label='Elastin')
    plt.plot(sv_stretch_var[64:n_zoom], sv_stress_var_collagen[64:n_zoom], '-.', linewidth=2, label='Collagen Total')
    plt.plot(sv_stretch_var[44:n_zoom], sv_stress_var_collagen_me[44:n_zoom], '+', linewidth=2, markersize=8, label='Collagen ME')
    plt.plot(sv_stretch_var[:n_zoom], sv_stress_var_collagen_ad[:n_zoom], '.', linewidth=2, markersize=8, label='Collagen AD')
    plt.plot(sv_stretch_var[44:n_zoom], sv_stress_var_muscle[44:n_zoom], 'x', linewidth=2, markersize=8, label='Muscle')
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
    else: 
        fig = ax.get_figure()

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
    ax.plot(time, results["muscle_cells"], color='black', linestyle='-', markevery=300, label='Muscle Cells',linewidth=1)

    if title:
        ax.set_title(title, fontsize=24, weight='bold')
    
    if xlabel:
        ax.set_xlabel('Time (years)', fontsize=20)
    if ylabel:
        ax.set_ylabel('Normalized Density', fontsize=20)
    
    ax.tick_params(axis='both', which='major', labelsize=18)

    ax.set_xlim(38, 75)
    ax.axvline(40, color='black', linestyle='--', linewidth=1.5)
    # ax.axvline(50, color='black', linestyle=':',  linewidth=1.5)
    ax.grid(True, linestyle='--', alpha=0.4)

    if legend:
        ax.legend(loc='upper right', fontsize=10, ncol=2)
    return fig

def plot_normalised_densities_by_genotype(results_tt, results_tc, results_cc):
    """
    Plot normalised densities for each genotype (TT, TC, CC) in a 2x2 grid.
    """
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])

    ax1 = fig.add_subplot(gs[0, 0])  # Top left: TT
    ax_legend = fig.add_subplot(gs[0, 1])  # Top right: Legend placeholder
    ax2 = fig.add_subplot(gs[1, 0])  # Bottom left: TC
    ax3 = fig.add_subplot(gs[1, 1])  # Bottom right: CC

    plot_normalised_densities(results_tt, ax1, title="Genotype: TT", ylabel=True, xlabel=False, legend=False)
    plot_normalised_densities(results_tc, ax2, title="Genotype: TC", ylabel=True, xlabel=True, legend=False)
    plot_normalised_densities(results_cc, ax3, title="Genotype: CC", ylabel=False, xlabel=True, legend=False)
    ax_legend.axis('off')
    handles, labels = ax1.get_legend_handles_labels()
    ax_legend.legend(handles, labels, loc='lower left', fontsize=10)

    plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.12, right=0.95, bottom=0.1, top=0.95)
    plt.show()

    return fig

def plot_normalised_densities_by_genotype2(results_tt, results_tc, results_cc):
    """
    Accessible side-by-side plot for genotypes with high-visibility markers and fonts.
    """
    fig, axes = plt.subplots(1, 3, figsize=(22, 8), sharey=True)
    
    # Plot each genotype
    plot_normalised_densities(results_tt, axes[0], title="Genotype: TT", ylabel=True, xlabel=True)
    plot_normalised_densities(results_tc, axes[1], title="Genotype: TC", ylabel=False, xlabel=True)
    plot_normalised_densities(results_cc, axes[2], title="Genotype: CC", ylabel=False, xlabel=True)
    
    # Collect legend handles from the first axis
    handles, labels = axes[0].get_legend_handles_labels()

    fig.legend(
        handles, labels,
        loc='lower center',
        fontsize=15,
        ncol=5,
        bbox_to_anchor=(0.5, -0.02),
        frameon=True,
        edgecolor='black',
        shadow=True
    )
    
    # axes[0].text(40.5, 1.6, "Infiltration", fontsize=16, rotation=90, verticalalignment='center', weight='bold')

    plt.show()
    
    return fig

def plot_normalised_densities_by_score(results_score_0, results_score_1, results_score_2, results_score_3, results_score_4):
    """
    Plot normalised densities for each polygenic score (0-4) in a 3x2 grid.
    """
    fig, axes = plt.subplots(3, 2, figsize=(18, 18))
    axes = axes.flatten()

    plot_normalised_densities(results_score_0, axes[0], title="Polygenic Score: 0", ylabel=True, xlabel=False, legend=False)
    plot_normalised_densities(results_score_1, axes[1], title="Polygenic Score: 1", ylabel=True, xlabel=False, legend=False)
    plot_normalised_densities(results_score_2, axes[2], title="Polygenic Score: 2", ylabel=True, xlabel=True, legend=False)
    plot_normalised_densities(results_score_3, axes[3], title="Polygenic Score: 3", ylabel=True, xlabel=True, legend=False)
    plot_normalised_densities(results_score_4, axes[4], title="Polygenic Score: 4", ylabel=True, xlabel=True, legend=False)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', fontsize=10)

    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()

    return fig

def plot_normalised_densities_by_score2(results_score_0, results_score_1, results_score_2, results_score_3, results_score_4):
    """
    Plot normalised densities for each polygenic score (0-4) side by side.
    """
    fig, axes = plt.subplots(1, 5, figsize=(20, 6))
    
    plot_normalised_densities(results_score_0, axes[0], title="Polygenic Score: 0", ylabel=True, xlabel=False, legend=False)
    plot_normalised_densities(results_score_1, axes[1], title="Polygenic Score: 1", ylabel=False, xlabel=True, legend=False)
    plot_normalised_densities(results_score_2, axes[2], title="Polygenic Score: 2", ylabel=False, xlabel=True, legend=False)
    plot_normalised_densities(results_score_3, axes[3], title="Polygenic Score: 3", ylabel=False, xlabel=True, legend=False)
    plot_normalised_densities(results_score_4, axes[4], title="Polygenic Score: 4", ylabel=False, xlabel=True, legend=False)
    
    for ax in axes:
        ax.set_ylim(-0.1, 1.8)
    
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', fontsize=10, ncol=7)

    plt.tight_layout(rect=[0, 0.08, 1, 1])  
    plt.show()
    
    return fig

def plot_stretch_by_genotype(results_tt, results_tc, results_cc):
    """
    Plot the systolic stretch over time for each genotype.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.plot(results_tt["time"], results_tt["lambda_sys"], label="TT", color='red')
    ax.plot(results_tc["time"], results_tc["lambda_sys"], label="TC", color='brown')
    ax.plot(results_cc["time"], results_cc["lambda_sys"], label="CC", color='gold')
    
    ax.set_xlabel("Time (years)", fontsize=20)
    ax.set_ylabel(r"Systolic Stretch $\lambda_{sys}$ (mm)", fontsize=20)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_xlim(40, 75)
    ax.legend(fontsize=12)
    ax.axvline(45, color='black', linestyle='--', linewidth=1.5)
    fig.tight_layout()
    plt.show()
    
    return fig

def plot_stretch_treat_notreat(results_tt, results_tc, results_cc, 
                               results_tt_treat, results_tc_treat, results_cc_treat):
    
    colors = ['#0072B2', '#D55E00', '#56B4E9'] 
    styles = ['-', '--', ':']
    
    fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(16, 10))
    
    genotypes = [(results_tt, results_tt_treat), (results_tc, results_tc_treat), (results_cc, results_cc_treat)]
    labels = ["TT", "TC", "CC"]

    for i, (res, res_treat) in enumerate(genotypes):
        ax[0].plot(res["time"], res["lambda_sys"], label=labels[i], 
                   color=colors[i], linestyle=styles[i], linewidth=4.5)
        ax[1].plot(res_treat["time"], res_treat["lambda_sys"], label=labels[i], 
                   color=colors[i], linestyle=styles[i], linewidth=4.5)

    for i, a in enumerate(ax):
        a.grid(True, linestyle='--', alpha=0.5, linewidth=1.5)
        a.tick_params(axis='both', labelsize=20)
        a.set_xlim(38, 75)
        
        # Consistent Y-axis padding for the bottom labels
        y_bottom = a.get_ylim()[0]
        y_label_pos = y_bottom + (a.get_ylim()[1] - y_bottom) * 0.05

        if i == 0:
            # Natural History Plot: Just Infiltration at x=40
            a.axvline(40, color='black', linestyle='-', linewidth=3, alpha=0.8)
            a.annotate("Infiltration Event", 
                       xy=(40, y_bottom), 
                       xytext=(41, y_label_pos), 
                       fontsize=18, 
                       weight='bold',
                       color='black',
                       bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="black", lw=2, alpha=1.0),
                       arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='black', lw=2.5))
        else: 
            # Treatment Plot: Both markers, label for Treatment at x=45
            a.axvline(40, color='black', linestyle='-', linewidth=3, alpha=0.4) # Faded to prioritize treatment line
            a.axvline(45, color='black', linestyle='-', linewidth=3, alpha=0.8)
            a.annotate("Treatment Applied",
                       xy=(45, y_bottom), 
                       xytext=(46, y_label_pos), 
                       fontsize=18, 
                       weight='bold',
                       color='black',
                       bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="black", lw=2, alpha=1.0),
                       arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='black', lw=2.5))

    ax[0].set_title('Without Intervention', fontsize=24, weight='bold', pad=20)
    ax[1].set_title('With Intervention', fontsize=24, weight='bold', pad=20)
    
    fig.text(0.5, 0.02, 'Time (Years)', ha='center', fontsize=26, weight='bold')
    fig.text(0.02, 0.5, r"Systolic Stretch $\lambda_{sys}$", va='center', rotation='vertical', fontsize=26, weight='bold')
    fig.suptitle('Influence of TGF-β Production on Arterial Stability', fontsize=30, weight='bold', y=0.98)
    
    ax[1].legend(fontsize=20, frameon=True, facecolor='white', edgecolor='black', loc='lower right', shadow=True, title='TGF-β Genotypes', title_fontsize=15)
    
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.93])
    plt.show()
    
    return fig

def plot_stretch_by_score(results_dict):
    """
    results_dict: {score_id: [list of simulation results]}
    Plots mean trajectory with shaded robust variance (95% CI) for polygenic scores.
    """
    fig, ax = plt.subplots(figsize=(16, 10))
    cmap = plt.get_cmap('viridis')
    
    # Sort keys to ensure color consistency (Score 0 to 4)
    for i, score in enumerate(sorted(results_dict.keys())):
        simulations = results_dict[score]
        color = cmap(i / 4.0)
        time = simulations[0]["time"]
        
        # Aggregate data and calculate robust statistics
        data_cube = np.array([sim["lambda_sys"] for sim in simulations])
        
        # Using Median and Percentiles to manage the non-linear "rupture" outliers
        center_line = np.median(data_cube, axis=0)
        lower_bound = np.percentile(data_cube, 2.5, axis=0)
        upper_bound = np.percentile(data_cube, 97.5, axis=0)
        
        # Plot variance (shaded area representing 95% Confidence Interval)
        ax.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.15)
        
        # Plot Trend Line (Median)
        ax.plot(time, center_line, label=f"Score {score}", 
                color=color, linewidth=5.0)

    # Styling and Grid
    ax.grid(True, linestyle='--', alpha=0.5, linewidth=1.5)
    ax.tick_params(axis='both', labelsize=20)
    ax.set_xlim(38, 75)
    
    # Vertical Line for Infiltration at t=40
    ax.axvline(40, color='black', linestyle='-', linewidth=2, alpha=0.8)

    # # --- ANNOTATION START ---
    # y_min, y_max = ax.get_ylim()
    # # Position the text box slightly above the bottom axis
    # y_label_pos = y_min + (y_max - y_min) * 0.05

    # ax.annotate("Infiltration Event", 
    #            xy=(40, y_min), 
    #            xytext=(38.5, y_label_pos), 
    #            fontsize=18, 
    #            weight='bold',
    #            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="black", lw=2),
    #            arrowprops=dict(arrowstyle="->", color='black', lw=2.5))
    # # --- ANNOTATION END ---
    
    # Typography
    ax.set_xlabel("Time (years)", fontsize=26, weight='bold', labelpad=15)
    ax.set_ylabel(r"Systolic Stretch $\lambda_{sys}$", fontsize=26, weight='bold', labelpad=15)
    ax.set_title('Influence of vSMC Volume on Arterial Stability', fontsize=30, weight='bold', pad=25)
    
    # Legend
    ax.legend(fontsize=16, loc='lower right', frameon=True, shadow=True, 
              edgecolor='black', title='Polygenic Scores', title_fontsize=16)
    
    fig.tight_layout()
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

def plot_elastin_smc(results_smc, results_no_smc):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        xlabel, ylabel = True, True

    time = results_smc["time"]

    ax.plot(time, results_smc["elastin_me"],   color='brown',   linestyle=':', label='Medial Collagen',    linewidth=1)
    ax.plot(time, results_no_smc["elastin_me"], linestyle='-', label='Medial Collagen',    linewidth=1)
    ax.legend(loc='upper right', fontsize=10, ncol=2)

    return ax

def plot_diameter_vs_time(results_tt, results_tc, results_cc):
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(results_tt["time"], results_tt["diameter"], label="TT", color='red')
    ax.plot(results_tc["time"], results_tc["diameter"], label="TC", color='brown')
    ax.plot(results_cc["time"], results_cc["diameter"], label="CC", color='gold')
    ax.set_xlabel("Time (years)", fontsize=20)
    ax.set_ylabel("Diameter (mm)", fontsize=20)   

    plt.show()
    return fig
            
def plot_stiffness_over_time(results_dict, params_dict):
    """
    results_dict: dict mapping score to results, e.g., {0: res0, 4: res4}
    params_dict: dict mapping score to corresponding ArterialParameters object
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    colors = {0: 'red', 1: 'brown', 2: 'gold', 3: 'green', 4: 'blue'}
    
    for score, results in results_dict.items():
        params = params_dict[score]
        time = results['time']
        stiffness = np.zeros_like(time)
        
        for i in range(len(time)):
            lam = results['lambda_sys'][i]
            
            # Perturb the stretch by 1% to measure the stiffness
            lam_pert = lam * 1.01 
            
            mE_M = results['elastin_me'][i]
            mC_M = results['collagen_me'][i]
            mC_A = results['collagen_ad'][i]
            mM = results['muscle_cells'][i]
            
            # Calculate the pressure difference using your force balance equation
            # force_balance_equation returns (calculated_pressure - target_pressure)
            res_base = functions.force_balance_equation([lam], mE_M, mC_M, mC_A, mM, params)
            res_pert = functions.force_balance_equation([lam_pert], mE_M, mC_M, mC_A, mM, params)
            
            # Change in pressure (dP) / Change in stretch (dlam)
            dP = res_pert - res_base
            dlam = lam_pert - lam
            
            # Store stiffness in kPa
            stiffness[i] = (dP / dlam) / 1000.0 
            
        ax.plot(time, stiffness, label=f"Score {score}", color=colors.get(score, 'black'), linewidth=2)

    ax.set_title('Arterial Stiffness (Tangent Modulus) Over Time', fontsize=16, weight='bold')
    ax.set_xlabel('Time (years)', fontsize=14)
    ax.set_ylabel('Stiffness / Tangent Modulus (kPa)', fontsize=14)
    ax.axvline(40, color='black', linestyle='--', linewidth=1.5, label='Immune Infiltration (t=40)')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(fontsize=12, loc='upper left')
    ax.set_xlim([35, 90])
    
    plt.show()
    return fig


def plot_sobol_indices(si_results, width=0.35):
    """
    Plot the Sobol sensitivity indices for the final circumferential stretch.
    """
    param_names = [
        'SMC Volume Fraction',
        'TGF-beta Level'
    ]
    s1 = si_results['S1']
    s2 = si_results['S2']
    st = si_results['ST']
    s1_error = si_results['S1_conf']
    s2_error = si_results['S2_conf']
    st_error = si_results['ST_conf']
    x = np.arange(len(param_names))
    print(x)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(
        x - width/2,
        s1,
        width,
        label='First-order (S1)',
        color='red'
    )
    ax.bar(
        x + width/2,
        st,
        width,
        label='Total-order (ST)',
        color='orange'
    )
    ax.errorbar(
        x - width/2,
        s1,
        yerr=s1_error,
        fmt='none',
        ecolor='black',
        capsize=5
    )
    ax.errorbar(
        x + width/2,
        st,
        yerr=st_error,
        fmt='none',
        ecolor='black',
        capsize=5
    )
    ax.set_title(
        'Sobol Sensitivity Indices for Final Circumferential Stretch',
        fontsize=14,
        weight='bold'
    )
    ax.set_xticks(x)
    ax.set_xticklabels(
        param_names,
        rotation=45,
        ha='right',
        fontsize=12
    )
    ax.set_ylabel('Sensitivity Index', fontsize=12)
    ax.grid(axis='y', linestyle='--', alpha=0.6)
    ax.legend(fontsize=11)

    plt.tight_layout()
    plt.show()

    return fig

def plot_stretch_heatmap(): 
    pass

def plot_genetic_risk_landscape(landscape_results):
    """
    Plots the pre-computed heatmap of the genetic risk landscape.
    Expects 'landscape_results' to contain: smc_range, tgf_range, and Z (2D array).
    """
    smc_range = landscape_results["smc_range"]
    tgf_range = landscape_results["tgf_range"]
    Z = landscape_results["Z"]

    norm = TwoSlopeNorm(vcenter=2.9, vmin=2.5, vmax=3.5) 
    
    plt.figure(figsize=(12, 9))
    im = plt.pcolormesh(smc_range, tgf_range, Z, cmap='RdBu_r', norm=norm, shading='gouraud')
    
    contours = plt.contour(smc_range, tgf_range, Z, levels=[2.88, 2.92], colors='white', alpha=0.5)
    plt.clabel(contours, inline=True, fontsize=10, fmt="Healthy Zone")

    plt.title("Arterial Health Landscape: TGF-B Level vs SMC Volume Fraction", fontsize=16, weight='bold')
    plt.xlabel("SMC Volume Fraction", fontsize=12)
    plt.ylabel("TGF-beta Expression Level", fontsize=12)
    plt.colorbar(im, label="Final Artery Diameter (mm)")
    
    plt.tight_layout()
    plt.show()

def plot_intervention_window_heatmap(landscape_data):
    """
    Plots a heatmap of final arterial diameter based on polygenic score and treatment timing.
    """
    genotype = landscape_data["genotype"]
    treat_times = landscape_data["treat_times"]
    final_diameters = landscape_data["final_diameters"]

    fig, ax = plt.subplots(figsize=(14, 8))
    vmax = max(np.max(final_diameters), 4.5)
    norm = TwoSlopeNorm(vcenter=2.9, vmin=2.5, vmax=vmax) # center around 2.9mm (healthy)

    im = ax.pcolormesh(treat_times, genotype, final_diameters, cmap='RdYlBu_r', norm=norm, shading='nearest')

    contours = ax.contour(treat_times, genotype, final_diameters, levels=[3.2, 4.0], colors=['white', 'black'], linewidths=2)
    fmt = {3.2: 'Stable (<3.2mm)', 5.0: 'Unstable (>5.0mm)'}
    ax.clabel(contours, inline=True, fontsize=12, fmt=fmt, weight='bold')

    ax.set_title("Intervention Timing vs. Final Artery Diameter", fontsize=22, weight='bold', pad=20)
    ax.set_xlabel("Time of Treatment (Years)", fontsize=16, weight='bold', labelpad=15)
    ax.set_ylabel("TGF-B Expression Level", fontsize=16, weight='bold', labelpad=15)
    
    ax.set_yticks(genotype)
    ax.tick_params(axis='both', labelsize=14)
    ax.axvline(x=40, color='red', linestyle='--', linewidth=3, alpha=0.8) # immune infiltration event
    ax.text(40.5, 3.8, 'Immune Infiltration\n(t=40)', color='red', fontsize=12, weight='bold', 
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='red', boxstyle='round,pad=0.3'))

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Final Artery Diameter at t=90 (mm)", fontsize=14, weight='bold', rotation=270, labelpad=25)
    cbar.ax.tick_params(labelsize=12)

    plt.tight_layout()
    plt.show()

    return fig

def plot_load_bearing_epochs(results_dict, params_dict, epochs=[35, 42, 55]):
    """
    Plots a stacked bar chart of load-bearing structures at specific time epochs.
    results_dict: dict mapping score to results {0: res0, 1: res1, ...}
    params_dict: dict mapping score to ArterialParameters
    """
    colors = {'SMC': '#D55E00', 'Elastin': '#E69F00', 'Collagen': '#0072B2'}
    
    fig, axes = plt.subplots(1, len(epochs), figsize=(18, 7), sharey=True)
    scores = sorted(list(results_dict.keys()))

    for ax, target_t in zip(axes, epochs):
        p_smc_list, p_el_list, p_col_me_list, p_col_ad_list = [], [], [], []

        for s in scores:
            res = results_dict[s]
            params = params_dict[s]
            
            # Find the time index closest to our target epoch
            idx = np.argmin(np.abs(res['time'] - target_t))
            
            # Extract variables at this specific time step
            lam = res['lambda_sys'][idx]
            mE = res['elastin_me'][idx]
            mC_m = res['collagen_me'][idx]
            mC_a = res['collagen_ad'][idx]
            mM = res['muscle_cells'][idx]
            
            # Calculate raw stresses using your physics engine
            sig_el = functions.sigma_elastin(lam, params)
            sig_smc = functions.sigma_muscle_t(lam, params)
            sig_col_m = functions.sigma_collagen_me(lam, params)
            sig_col_a = functions.sigma_collagen_ad(lam, params)
            
            # Convert Cauchy stress to transmural pressure components (kPa)
            prefactor = 1 / (params.c_radius_tzero * lam**2 * params.c_lambda_z)
            
            load_el = prefactor * params.c_thickness_me * mE * sig_el
            load_smc = prefactor * params.c_thickness_me * mM * sig_smc
            load_col_me = prefactor * (params.c_thickness_me * mC_m * sig_col_m)
                                    
            load_col_ad = prefactor * (params.c_thickness_ad * mC_a * sig_col_a)
            
            p_el_list.append(load_el / 1000)
            p_smc_list.append(load_smc / 1000)
            p_col_me_list.append(load_col_me / 1000)
            p_col_ad_list.append(load_col_ad / 1000)
            
        # Create Stacked Bars
        x = np.arange(len(scores))
        width = 0.6
        
        p_el, p_smc, p_col_me, p_col_ad = np.array(p_el_list), np.array(p_smc_list), np.array(p_col_me_list), np.array(p_col_ad_list)
        
        ax.bar(x, p_el, width, label='Elastin', color=colors['Elastin'], edgecolor='black', linewidth=1.2)
        ax.bar(x, p_smc, width, bottom=p_el, label='VSMCs', color=colors['SMC'], edgecolor='black', linewidth=1.2)
        ax.bar(x, p_col_me, width, bottom=p_el + p_smc, label='Collagen (Media)', color=colors['Collagen'], edgecolor='black', linewidth=1.2)
        ax.bar(x, p_col_ad, width, bottom=p_el + p_smc + p_col_me, label='Collagen (Adventitia)', color='#009E73', edgecolor='black', linewidth=1.2)

        
        # Poster-ready typography
        epoch_labels = {35: "Pre-Infiltration\n(Homeostasis, t=35)", 
                        42: "Active Remodeling\n(Post-Infiltration, t=42)", 
                        55: "Late Stage\n(Degeneration, t=55)"}
        
        ax.set_title(epoch_labels.get(target_t, f"Year {target_t}"), fontsize=18, weight='bold', pad=15)
        ax.set_xticks(x)
        ax.set_xticklabels([f"Score {s}" for s in scores], fontsize=16, weight='bold')
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        ax.tick_params(axis='y', labelsize=14)

    # Global labels and legend
    axes[0].set_ylabel("Load Borne (kPa)", fontsize=20, weight='bold')
    fig.suptitle('Artery Wall Load Transfer Mechanics by Polygenic Score', fontsize=26, weight='bold', y=0.98)
    
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=18, bbox_to_anchor=(0.5, -0.02), frameon=True, edgecolor='black')
    
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    plt.show()
    return fig