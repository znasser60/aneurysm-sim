import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.transforms as mtransforms
import matplotlib.patheffects as path_effects

from aneurysm_sim.model import model


def plot_pressure_vs_stretch(results, n_zoom=120):
    """
    Plots the pressure vs stretch of the artery with a zoomed-in view.

    Parameters:
    results (dict):
        A dictionary containing the stress/pressure simulation results.
    n_zoom (int):
        The number of data points to display in the zoomed-in view.
    """
    stretch = results["sv_stretch_var"]
    sv_pressure_var = results["sv_pressure_var"]
    sv_pressure_var_elastin = results["sv_pressure_var_elastin"]
    sv_pressure_var_collagen = results["sv_pressure_var_collagen"]
    sv_pressure_var_muscle_p = results["sv_pressure_var_muscle_p"]
    sv_pressure_var_muscle_a = results["sv_pressure_var_muscle_a"]

    plt.figure(figsize=(12, 8))
    plt.plot(
        stretch[:n_zoom], sv_pressure_var[:n_zoom] / 1e3, linewidth=2, label="Total"
    )
    plt.plot(
        stretch[32:n_zoom],
        sv_pressure_var_elastin[32:n_zoom] / 1e3,
        "--",
        linewidth=2,
        label="Elastin",
    )
    plt.plot(
        stretch[64:n_zoom],
        sv_pressure_var_collagen[64:n_zoom] / 1e3,
        "-.",
        linewidth=2,
        label="Collagen",
    )
    plt.plot(
        stretch[44:n_zoom],
        sv_pressure_var_muscle_p[44:n_zoom] / 1e3,
        "--",
        linewidth=2,
        label="Muscle Passive",
    )
    plt.plot(
        stretch[:n_zoom],
        sv_pressure_var_muscle_a[:n_zoom] / 1e3,
        "--",
        linewidth=2,
        label="Muscle Active",
    )
    plt.title("Pressure vs Stretch")
    plt.xlabel("Stretch")
    plt.ylabel("Pressure (kPa)")
    plt.legend(loc="upper left")
    plt.grid(True)
    plt.ylim([0, 60])
    plt.show()


def plot_pressure_vs_diameter(results, n_zoom=120):
    """
    Plots the pressure vs diameter of the artery with a zoomed-in view.

    Parameters:
    results (dict):
        A dictionary containing the stress/pressure simulation results.
    n_zoom (int): 
        The number of data points to display in the zoomed-in view.
    """
    sv_diam_var = results["sv_diam_var"]
    sv_pressure_var = results["sv_pressure_var"]
    sv_pressure_var_elastin = results["sv_pressure_var_elastin"]
    sv_pressure_var_collagen = results["sv_pressure_var_collagen"]
    sv_pressure_var_collagen_me = results["sv_pressure_var_collagen_me"]
    sv_pressure_var_collagen_ad = results["sv_pressure_var_collagen_ad"]
    sv_pressure_var_muscle_p = results["sv_pressure_var_muscle_p"]
    sv_pressure_var_muscle_a = results["sv_pressure_var_muscle_a"]
    sv_pressure_var_muscle = results["sv_pressure_var_muscle"]

    plt.figure(figsize=(12, 8))
    plt.plot(
        sv_diam_var[:n_zoom],
        sv_pressure_var[:n_zoom] / 1e3,
        "-",
        linewidth=6,
        label="Total",
    )
    plt.plot(
        sv_diam_var[32:n_zoom],
        sv_pressure_var_elastin[32:n_zoom] / 1e3,
        "--",
        linewidth=3,
        label="Elastin",
    )
    plt.plot(
        sv_diam_var[64:n_zoom],
        sv_pressure_var_collagen[64:n_zoom] / 1e3,
        "-.",
        linewidth=3,
        label="Collagen Total",
    )
    plt.plot(
        sv_diam_var[64:n_zoom],
        sv_pressure_var_collagen_me[64:n_zoom] / 1e3,
        "-.",
        linewidth=3,
        label="Collagen ME",
    )
    plt.plot(
        sv_diam_var[64:n_zoom],
        sv_pressure_var_collagen_ad[64:n_zoom] / 1e3,
        "-.",
        linewidth=3,
        label="Collagen AD",
    )
    plt.plot(
        sv_diam_var[44:n_zoom],
        sv_pressure_var_muscle[44:n_zoom] / 1e3,
        "+",
        linewidth=2,
        markersize=8,
        label="vSMC Total",
    )
    plt.plot(
        sv_diam_var[44:n_zoom],
        sv_pressure_var_muscle_p[44:n_zoom] / 1e3,
        "+",
        linewidth=2,
        markersize=8,
        label="vSMCp",
    )
    plt.plot(
        sv_diam_var[:n_zoom],
        sv_pressure_var_muscle_a[:n_zoom] / 1e3,
        "+",
        linewidth=2,
        markersize=8,
        label="vSMCa",
    )
    plt.axvline(x=2.9, color="red", linestyle="-", linewidth=2)
    plt.axhline(y=16, color="red", linestyle="-", linewidth=2)
    plt.plot(2.9, 16, "ro", linewidth=3, markersize=8)

    # plt.title('Pressure vs Diameter', fontsize=30, weight='bold', pad=25)
    plt.xlabel("Diameter (mm)", fontsize=20, weight="bold", labelpad=15)
    plt.ylabel("Pressure (kPa)", fontsize=20, weight="bold", labelpad=15)
    plt.legend(
        loc="upper left", fontsize=16, frameon=True, shadow=True, edgecolor="black"
    )
    plt.ylim([0, 60])
    plt.tick_params(axis="both", which="major", labelsize=18)
    plt.grid(True)
    plt.show()


def plot_stretch_vs_stress(results, n_zoom=120):
    """
    Plots the stretch vs stress with a zoomed-in view.

    Parameters:
    results (dict):
        A dictionary containing the stress/pressure simulation results.
    n_zoom (int):
        The number of data points to display in the zoomed-in view.
    """
    sv_stretch_var = results["sv_stretch_var"]
    sv_stress_var_elastin = results["sv_stress_var_elastin"]
    sv_stress_var_collagen_me = results["sv_stress_var_collagen_me"]
    sv_stress_var_collagen_ad = results["sv_stress_var_collagen_ad"]
    sv_stress_var_collagen = results["sv_stress_var_collagen"]
    sv_stress_var_muscle = results["sv_stress_var_muscle_t"]
    sv_stress_var_muscle_a = results["sv_stress_var_muscle_a"]
    sv_stress_var_muscle_p = results["sv_stress_var_muscle_p"]
    sv_stress_var_total = results["sv_stress_var_total"]

    plt.figure(figsize=(12, 8))
    ax = plt.gca()
    plt.plot(
        sv_stretch_var[32:n_zoom],
        sv_stress_var_elastin[32:n_zoom],
        "--",
        linewidth=2,
        label="Elastin",
    )
    plt.plot(
        sv_stretch_var[64:n_zoom],
        sv_stress_var_collagen[64:n_zoom],
        "-.",
        linewidth=2,
        label="Collagen",
    )
    plt.plot(
        sv_stretch_var[44:n_zoom],
        sv_stress_var_collagen_me[44:n_zoom],
        "+",
        linewidth=2,
        markersize=8,
        label="Collagen ME",
    )
    plt.plot(
        sv_stretch_var[:n_zoom],
        sv_stress_var_collagen_ad[:n_zoom],
        ".",
        linewidth=2,
        markersize=8,
        label="Collagen AD",
    )
    plt.plot(
        sv_stretch_var[64:n_zoom],
        sv_stress_var_muscle[64:n_zoom],
        "-",
        linewidth=2,
        markersize=8,
        label="vSMC",
    )
    plt.plot(
        sv_stretch_var[44:n_zoom],
        sv_stress_var_muscle_a[44:n_zoom],
        "+",
        linewidth=2,
        markersize=3,
        label="aVSMC",
    )
    plt.plot(
        sv_stretch_var[44:n_zoom],
        sv_stress_var_muscle_p[44:n_zoom],
        ".",
        linewidth=2,
        markersize=3,
        label="pVSMC",
    )

    plt.plot(
        sv_stretch_var[:n_zoom],
        sv_stress_var_total[:n_zoom],
        "-",
        linewidth=3,
        label="Total Stress",
    )
    # plt.title('Component-specific Stretch vs. Stress', fontsize=30, weight='bold', pad=25)
    plt.axvspan(1.0, 1.3, alpha=0.12, color="gray", zorder=0)
    plt.axvspan(1.3, 1.8, alpha=0.12, color="red", zorder=0)
    plt.xlim(0.6, 1.7)
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(
        1.15,
        0.95,
        "Normal Physiological",
        transform=trans,
        ha="center",
        va="top",
        fontsize=14,
        color="gray",
        weight="bold",
    )
    ax.text(
        1.5,
        0.95,
        "Pathological",
        transform=trans,
        ha="center",
        va="top",
        fontsize=14,
        color="firebrick",
        weight="bold",
    )
    plt.xlabel(r"Stretch $\lambda$", fontsize=20, weight="bold", labelpad=15)
    plt.ylabel(r"Stress $\sigma$", fontsize=20, weight="bold", labelpad=15)
    plt.tick_params(axis="both", which="major", labelsize=18)
    plt.legend(
        loc="upper left", fontsize=16, frameon=True, shadow=True, edgecolor="black"
    )
    plt.grid(True)
    plt.show()


def plot_normalised_densities(
    results, ax=None, title=None, legend=False, xlabel=False, ylabel=False
):
    """
    Plots the normalised densities of elastin, collagen, immune cells, latent and active TGF Beta.

    Parameters:
    results (dict): 
        A dictionary containing the aneurysm simultion results.
    Remaining parameters are for customizing the plot appearance.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        xlabel, ylabel = True, True
    else:
        fig = ax.get_figure()

    time = results["time"]

    ax.plot(
        time,
        results["collagen_me"],
        color="brown",
        linestyle=":",
        label="Medial Collagen",
        linewidth=1,
    )
    ax.plot(
        time,
        results["elastin_me"],
        color="brown",
        linestyle=":",
        label="Medial Elastin",
        linewidth=1,
    )
    ax.plot(
        time,
        results["elastases"],
        color="red",
        marker="o",
        markevery=300,
        label="Elastases",
        linewidth=1,
    )
    ax.plot(
        time,
        results["collagenases"],
        color="red",
        marker="o",
        markevery=300,
        label="Collagenases",
        linewidth=1,
    )
    ax.plot(
        time, results["fibroblast"], color="brown", linewidth=3, label="Fibroblasts"
    )
    ax.plot(
        time,
        results["collagen_ad"],
        color="brown",
        linewidth=1,
        label="Adventitial Collagen",
    )
    ax.plot(
        time, results["procollagen"], color="brown", linewidth=3, label="Procollagen"
    )
    ax.plot(
        time,
        results["collagenase"],
        color="magenta",
        marker="^",
        markevery=300,
        linestyle="-",
        label="Collagenase",
        linewidth=1,
    )
    ax.plot(
        time,
        results["zymogen"],
        color="magenta",
        marker="^",
        markevery=300,
        linestyle="--",
        label="Zymogen",
        linewidth=1,
    )
    ax.plot(
        time,
        results["timp"],
        color="gold",
        marker="v",
        markevery=300,
        linestyle="-",
        label="TIMP",
        linewidth=1,
    )
    ax.plot(
        time,
        results["immune_cells"],
        color="red",
        marker="o",
        markevery=300,
        label="Immune Cells",
        linewidth=3,
    )
    ax.plot(
        time,
        results["latent_tgf_beta"],
        color="green",
        linestyle="--",
        marker="s",
        markevery=300,
        label="Latent TGF-β",
        linewidth=1,
    )
    ax.plot(
        time,
        results["active_tgf_beta"],
        color="green",
        linestyle="-",
        marker="s",
        markevery=300,
        label="Active TGF-β",
        linewidth=1,
    )
    ax.plot(
        time,
        results["muscle_cells"],
        color="black",
        linestyle="-",
        markevery=300,
        label="Muscle Cells",
        linewidth=1,
    )

    if title:
        ax.set_title(title, fontsize=24, weight="bold")

    if xlabel:
        ax.set_xlabel("Time (years)", fontsize=20, weight="bold", labelpad=15)
    if ylabel:
        ax.set_ylabel("Normalized Density", fontsize=20, weight="bold", labelpad=15)

    ax.tick_params(axis="both", which="major", labelsize=18)

    ax.set_xlim(38, 75)
    ax.axvline(40, color="black", linestyle="--", linewidth=1.5)
    ax.grid(True, linestyle="--", alpha=0.4)

    if legend:
        ax.legend(loc="upper right", fontsize=10, ncol=2)
    return fig


def plot_normalised_densities_by_genotype(results_tt, results_tc, results_cc):
    """
    Plots the normalised densities of elastin, collagen, immune cells, latent 
    and active TGF Beta for each genotype side by side.

    Parameters:
    results_tt (dict): 
        A dictionary containing the aneurysm simulation results for genotype TT.
    results_tc (dict):
        A dictionary containing the aneurysm simulation results for genotype TC.
    results_cc (dict):
        A dictionary containing the aneurysm simulation results for genotype CC.
    """
    fig, axes = plt.subplots(1, 3, figsize=(22, 8), sharey=True)

    plot_normalised_densities(
        results_tt, axes[0], title="Genotype TT (0.7)", ylabel=True, xlabel=False
    )
    plot_normalised_densities(
        results_tc, axes[1], title="Genotype TC (0.9)", ylabel=False, xlabel=True
    )
    plot_normalised_densities(
        results_cc, axes[2], title="Genotype CC (1.1)", ylabel=False, xlabel=False
    )

    handles, labels = axes[0].get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        loc="lower center",
        fontsize=15,
        ncol=5,
        bbox_to_anchor=(0.5, -0.02),
        frameon=True,
        edgecolor="black",
        shadow=True,
    )

    plt.show()

    return fig


def plot_normalised_densities_by_score(
    results_score_0, results_score_1, results_score_2, results_score_3, results_score_4
):
    """
    Plot normalised densities for each polygenic score (0-4) in a 3x2 grid layout.

    Parameters:
    results_score_0 (dict): 
        A dictionary containing the aneurysm simulation results for polygenic score 0.
    results_score_1 (dict):
        A dictionary containing the aneurysm simulation results for polygenic score 1.
    results_score_2 (dict): 
        A dictionary containing the aneurysm simulation results for polygenic score 2.
    results_score_3 (dict):
        A dictionary containing the aneurysm simulation results for polygenic score 3.
    results_score_4 (dict):
        A dictionary containing the aneurysm simulation results for polygenic score 4.
    """
    fig, axes = plt.subplots(3, 2, figsize=(18, 18))
    axes = axes.flatten()

    plot_normalised_densities(
        results_score_0,
        axes[0],
        title="Polygenic Score: 0",
        ylabel=True,
        xlabel=False,
        legend=False,
    )
    plot_normalised_densities(
        results_score_1,
        axes[1],
        title="Polygenic Score: 1",
        ylabel=True,
        xlabel=False,
        legend=False,
    )
    plot_normalised_densities(
        results_score_2,
        axes[2],
        title="Polygenic Score: 2",
        ylabel=True,
        xlabel=True,
        legend=False,
    )
    plot_normalised_densities(
        results_score_3,
        axes[3],
        title="Polygenic Score: 3",
        ylabel=True,
        xlabel=True,
        legend=False,
    )
    plot_normalised_densities(
        results_score_4,
        axes[4],
        title="Polygenic Score: 4",
        ylabel=True,
        xlabel=True,
        legend=False,
    )

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", fontsize=10)

    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()

    return fig


def plot_normalised_densities_by_score_reduced(
    results_score_0, results_score_2, results_score_4
):
    """
    Plot normalised densities for polygenic scores 0, 2, and 4 in a single row layout.

    Parameters:
    results_score_0 (dict): 
        A dictionary containing the aneurysm simulation results for polygenic score 0.
    results_score_2 (dict):
        A dictionary containing the aneurysm simulation results for polygenic score 2.
    results_score_4 (dict):
        A dictionary containing the aneurysm simulation results for polygenic score 4.
    
    Note: This plot is for reduced polygenic scores to focus on the extremes. 
    """
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)

    plot_normalised_densities(
        results_score_0,
        axes[0],
        title="Score 0 (73%)",
        ylabel=True,
        xlabel=False,
        legend=False,
    )
    plot_normalised_densities(
        results_score_2,
        axes[1],
        title="Score 2 (67%)",
        ylabel=False,
        xlabel=True,
        legend=False,
    )
    plot_normalised_densities(
        results_score_4,
        axes[2],
        title="Score 4 (49%)",
        ylabel=False,
        xlabel=False,
        legend=False,
    )

    for ax in axes:
        ax.set_ylim(-0.1, 1.8)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        fontsize=15,
        ncol=5,
        bbox_to_anchor=(0.5, -0.02),
        frameon=True,
        edgecolor="black",
        shadow=True,
    )

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.show()

    return fig

def plot_stretch(results):
    """
    Plot patient-specific systolic stretch over time based on genotype input.

    Parameters:
    results (dict): 
        A dictionary containing the patient aneurysm simulation results.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    ax.plot(results["time"], results["lambda_sys"], color="blue")

    ax.set_xlabel("Time (years)", fontsize=20)
    ax.set_ylabel(r"Systolic Stretch $\lambda_{sys}$ (mm)", fontsize=20)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.set_xlim(40, 75)
    ax.legend(fontsize=12)
    ax.axvline(45, color="black", linestyle="--", linewidth=1.5)
    fig.tight_layout()
    plt.show()

    return fig

def plot_stretch_by_genotype(results_tt, results_tc, results_cc):
    """
    Plot the systolic stretch over time for each genotype.

    Parameters:
    results_tt (dict): 
        A dictionary containing the aneurysm simulation results for genotype TT.
    results_tc (dict):
        A dictionary containing the aneurysm simulation results for genotype TC.
    results_cc (dict):
        A dictionary containing the aneurysm simulation results for genotype CC.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    ax.plot(results_tt["time"], results_tt["lambda_sys"], label="TT", color="red")
    ax.plot(results_tc["time"], results_tc["lambda_sys"], label="TC", color="brown")
    ax.plot(results_cc["time"], results_cc["lambda_sys"], label="CC", color="gold")

    ax.set_xlabel("Time (years)", fontsize=20)
    ax.set_ylabel(r"Systolic Stretch $\lambda_{sys}$ (mm)", fontsize=20)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.set_xlim(40, 75)
    ax.legend(fontsize=12)
    ax.axvline(45, color="black", linestyle="--", linewidth=1.5)
    fig.tight_layout()
    plt.show()

    return fig


def plot_stretch_by_geno_treat(
    results_tt,
    results_tc,
    results_cc,
    results_tt_treat,
    results_tc_treat,
    results_cc_treat,
):
    """
    Plot the systolic stretch over time for each genotype, 
    comparing untreated and treated cases.

    Parameters:
    results_tt/results_tc/results_cc (dict): 
        A dictionary containing the aneurysm simulation results for each genotype (untreated).
    results_tt_treat/results_tc_treat/results_cc_treat (dict):
        A dictionary containing the aneurysm simulation results for each genotype (treated).
    """
    colors = ["#0072B2", "#D55E00", "#56B4E9"]
    styles = ["-", "--", ":"]

    fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(16, 10))

    genotypes = [
        (results_tt, results_tt_treat),
        (results_tc, results_tc_treat),
        (results_cc, results_cc_treat),
    ]
    labels = ["TT", "TC", "CC"]

    for i, (res, res_treat) in enumerate(genotypes):
        ax[0].plot(
            res["time"],
            res["lambda_sys"],
            label=labels[i],
            color=colors[i],
            linestyle=styles[i],
            linewidth=4.5,
        )
        ax[1].plot(
            res_treat["time"],
            res_treat["lambda_sys"],
            label=labels[i],
            color=colors[i],
            linestyle=styles[i],
            linewidth=4.5,
        )

    for i, a in enumerate(ax):
        a.grid(True, linestyle="--", alpha=0.5, linewidth=1.5)
        a.tick_params(axis="both", labelsize=20)
        a.set_xlim(38, 75)
        y_bottom = a.get_ylim()[0]
        y_label_pos = y_bottom + (a.get_ylim()[1] - y_bottom) * 0.05

        if i == 0:
            # Natural History Plot: Just Infiltration at x=40
            a.axvline(40, color="black", linestyle="-", linewidth=3, alpha=0.8)
            a.annotate(
                "Infiltration Event",
                xy=(40, y_bottom),
                xytext=(41, y_label_pos),
                fontsize=18,
                weight="bold",
                color="black",
                bbox=dict(
                    boxstyle="round,pad=0.4", fc="white", ec="black", lw=2, alpha=1.0
                ),
                arrowprops=dict(
                    arrowstyle="->", connectionstyle="arc3", color="black", lw=2.5
                ),
            )
        else:
            a.axvline(40, color="black", linestyle="-", linewidth=3, alpha=0.4)
            a.axvline(45, color="black", linestyle="-", linewidth=3, alpha=0.8)
            a.annotate(
                "Treatment Applied",
                xy=(45, y_bottom),
                xytext=(46, y_label_pos),
                fontsize=18,
                weight="bold",
                color="black",
                bbox=dict(
                    boxstyle="round,pad=0.4", fc="white", ec="black", lw=2, alpha=1.0
                ),
                arrowprops=dict(
                    arrowstyle="->", connectionstyle="arc3", color="black", lw=2.5
                ),
            )

    ax[0].set_title("Without Intervention", fontsize=24, weight="bold", pad=20)
    ax[1].set_title("With Intervention", fontsize=24, weight="bold", pad=20)

    fig.text(0.5, 0.02, "Time (Years)", ha="center", fontsize=26, weight="bold")
    fig.text(
        0.02,
        0.5,
        r"Systolic Stretch $\lambda_{sys}$",
        va="center",
        rotation="vertical",
        fontsize=26,
        weight="bold",
    )
    # fig.suptitle('Influence of TGF-β Production on Arterial Stability', fontsize=30, weight='bold', y=0.98)

    ax[1].legend(
        fontsize=20,
        frameon=True,
        facecolor="white",
        edgecolor="black",
        loc="lower right",
        shadow=True,
    )

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.93])
    plt.show()

    return fig


def plot_stretch_by_score(results_dict, results_dict_treat):
    """
    Plots the systolic stretch over time for each polygenic score, 
    comparing untreated and treated cases.

    Parameters:
    results_dict (dict):
        A dictionary containing the aneurysm simulation results for each polygenic score (untreated).
    results_dict_treat (dict):
        A dictionary containing the aneurysm simulation results for each polygenic score (treated).
    """
    fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(18, 10))

    cmap = plt.get_cmap("viridis")
    labels = {0: "Score 0", 1: "Score 1", 2: "Score 2", 3: "Score 3", 4: "Score 4"}

    panels = [(ax[0], results_dict), (ax[1], results_dict_treat)]

    for panel_ax, res_dict in panels:
        for i, score in enumerate(sorted(res_dict.keys())):
            simulations, mean_result = res_dict[score]
            color = cmap(i / 4.0)
            time = simulations[0]["time"]

            data_cube = np.array([sim["lambda_sys"] for sim in simulations])
            center_line = np.array(mean_result["lambda_sys"])
            lower_bound = np.min(data_cube, axis=0)
            upper_bound = np.max(data_cube, axis=0)

            panel_ax.fill_between(
                time, lower_bound, upper_bound, color=color, alpha=0.15
            )
            panel_ax.plot(
                time, center_line, label=labels[score], color=color, linewidth=5.0
            )

    for i, a in enumerate(ax):
        a.grid(True, linestyle="--", alpha=0.5, linewidth=1.5)
        a.tick_params(axis="both", labelsize=20)
        a.set_xlim(38, 75)
        y_bottom = a.get_ylim()[0]
        y_label_pos = y_bottom + (a.get_ylim()[1] - y_bottom) * 0.05

        if i == 0:
            # Natural History Plot: Just Infiltration at x=40
            a.axvline(40, color="black", linestyle="-", linewidth=3, alpha=0.8)
            a.annotate(
                "Infiltration Event",
                xy=(40, y_bottom),
                xytext=(41, y_label_pos),
                fontsize=18,
                weight="bold",
                color="black",
                bbox=dict(
                    boxstyle="round,pad=0.4", fc="white", ec="black", lw=2, alpha=1.0
                ),
                arrowprops=dict(
                    arrowstyle="->", connectionstyle="arc3", color="black", lw=2.5
                ),
            )
        else:
            a.axvline(40, color="black", linestyle="-", linewidth=3, alpha=0.4)
            a.axvline(45, color="black", linestyle="-", linewidth=3, alpha=0.8)
            a.annotate(
                "Treatment Applied",
                xy=(45, y_bottom),
                xytext=(46, y_label_pos),
                fontsize=18,
                weight="bold",
                color="black",
                bbox=dict(
                    boxstyle="round,pad=0.4", fc="white", ec="black", lw=2, alpha=1.0
                ),
                arrowprops=dict(
                    arrowstyle="->", connectionstyle="arc3", color="black", lw=2.5
                ),
            )

    ax[0].set_title("Without Intervention", fontsize=24, weight="bold", pad=20)
    ax[1].set_title("With Intervention", fontsize=24, weight="bold", pad=20)

    fig.text(0.5, 0.02, "Time (Years)", ha="center", fontsize=26, weight="bold")
    fig.text(
        0.02,
        0.5,
        r"Systolic Stretch $\lambda_{sys}$",
        va="center",
        rotation="vertical",
        fontsize=26,
        weight="bold",
    )

    ax[1].legend(
        fontsize=18,
        frameon=True,
        facecolor="white",
        edgecolor="black",
        loc="lower right",
        shadow=True,
        title="Polygenic Scores",
        title_fontsize=18,
    )

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.93])
    plt.show()

    return fig


def plot_stretch_evolution_dist(
    results,
    stretch_type="attachment",
    target_years=[40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 60],
):
    """
    Plots multiple triangular distributions layered on the same axes.
    Specifically targets requested years to capture active remodeling.

    Parameters:
    results (dict):
        A dictionary containing the aneurysm simulation results.
    stretch_type (str):
        Type of stretch to plot. Options are "attachment" or "recruitment".
    target_years (list):
        List of years to plot the distributions for. Default is 40-60. 
    """
    time = results["time"]

    if stretch_type == "attachment":
        l_min_raw = results["lambda_att_min"]
        l_mode_raw = results["lambda_att_mode"]
        l_max_raw = results["lambda_att_max"]
        xlabel = r"Attachment Stretch $\lambda_A^{AT}$"
    elif stretch_type == "recruitment":
        l_min_raw = results["lambda_rec_min"]
        l_mode_raw = results["lambda_rec_mode"]
        l_max_raw = results["lambda_rec_max"]
        xlabel = r"Recruitment Stretch $\lambda_A^R$"
    else:
        raise ValueError("stretch_type must be 'attachment' or 'recruitment'")

    indices = [np.argmin(np.abs(time - y)) for y in target_years]

    fig, ax = plt.subplots(figsize=(10, 6))

    cmap = plt.cm.viridis
    norm = mpl.colors.Normalize(vmin=min(target_years), vmax=max(target_years))

    all_selected_vals = np.concatenate(
        [l_min_raw[indices], l_mode_raw[indices], l_max_raw[indices]]
    )
    x_min_global = np.min(all_selected_vals) * 0.98
    x_max_global = np.max(all_selected_vals) * 1.02
    x = np.linspace(x_min_global, x_max_global, 1500)

    peaks_x = []
    peaks_y = []

    for i, idx in enumerate(indices):
        t = time[idx]
        vertices = np.sort([l_min_raw[idx], l_mode_raw[idx], l_max_raw[idx]])
        a, b, c = vertices
        if c - a < 1e-6:
            y = np.zeros_like(x)
            peak_y = 0
        else:
            y = np.where(
                (x >= a) & (x <= b),
                2 * (x - a) / ((b - a + 1e-9) * (c - a)),
                np.where(
                    (x > b) & (x <= c), 2 * (c - x) / ((c - b + 1e-9) * (c - a)), 0
                ),
            )
            peak_y = 2 / (c - a)

        color = cmap(norm(t))
        ax.fill_between(x, y, alpha=0.15, color=color, edgecolor="none")
        ax.plot(x, y, color=color, linewidth=2.5)
        peaks_x.append(b)
        peaks_y.append(peak_y)

    ax.set_xlabel(xlabel, fontsize=18, weight="bold", labelpad=15)
    ax.set_ylabel("Probability Density", fontsize=18, weight="bold", labelpad=15)
    ax.tick_params(axis="both", which="major", labelsize=14)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)

    ax.grid(True, linestyle="--", alpha=0.4)

    if len(peaks_x) >= 2:
        x_start, y_start = peaks_x[0], peaks_y[0]
        x_end, y_end = peaks_x[-1], peaks_y[-1]

        ax.annotate(
            "Time",
            xy=(x_end, y_end),
            xycoords="data",
            xytext=(x_start, y_start),
            textcoords="data",
            arrowprops=dict(
                arrowstyle="->",  
                color="black",
                lw=2.5,
                connectionstyle="arc3,rad=0.15",  
            ),
            fontsize=16,
            weight="bold",
            color="black",
            ha="center",
            va="bottom",
        )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, aspect=30, pad=0.03)
    cbar.set_label("Simulation Time (Years)", fontsize=16, weight="bold", labelpad=15)
    cbar.ax.tick_params(labelsize=14)

    plt.tight_layout()
    plt.show()


def plot_load_bearing(params_dict):
    """
    Grid of pie charts showing the design-target load-bearing partition.

    Parameters:
    params_dict (dict):
        A dictionary where keys are polygenic scores and values are calculated parameters. 
    """
    labels = ["Elastin", "VSMC", "Collagen (Media)", "Collagen (Adventitia)"]
    palette = {
        "Elastin": "#64DFDF", 
        "VSMC": "#6930C3", 
        "Collagen (Media)": "#5390D9", 
        "Collagen (Adventitia)": "#48BFE3", 
    }
    wedge_colors = [palette[l] for l in labels]
    explode = [0.0, 0.10, 0.0, 0.0] 

    scores = sorted(params_dict.keys())
    n = len(scores)

    plt.rcParams["font.family"] = "sans-serif"
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 4.8))
    if n == 1:
        axes = [axes]

    ink = "#1a1a2e"

    for ax, s in zip(axes, scores):
        p = params_dict[s]

        smc = p.c_load_borne_muscle_p + p.c_load_borne_muscle_a
        elastin = p.c_load_borne_elastin
        ratio = p.c_collagen_ratio_ad_me
        col_me = p.c_load_borne_collagen * (1 / (1 + ratio))
        col_ad = p.c_load_borne_collagen * (ratio / (1 + ratio))
        fractions = [elastin, smc, col_me, col_ad]

        wedges, texts, autotexts = ax.pie(
            fractions,
            explode=explode,
            colors=wedge_colors,
            autopct=lambda pct: f"{pct:.0f}%" if pct > 4 else "",
            pctdistance=0.72,
            startangle=90,
            counterclock=False,
            wedgeprops={"edgecolor": "white", "linewidth": 2.0, "antialiased": True},
            textprops={"color": "white", "fontsize": 12, "weight": "bold"},
        )

        wedges[1].set_path_effects(
            [path_effects.withSimplePatchShadow(offset=(2, -2), alpha=0.22)]
        )

        for at in autotexts:
            at.set_path_effects(
                [path_effects.withStroke(linewidth=2.5, foreground="#00000055")]
            )

        ax.set_title(f"Score {s}", fontsize=15, weight="bold", pad=12, color=ink)
        ax.text(
            0.5,
            -0.06,
            f"vSMC frac = {p.smc_fraction:.2f}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=10,
            color="#666",
            style="italic",
        )
        ax.set_aspect("equal")

    handles = [
        plt.Rectangle(
            (0, 0), 1, 1, facecolor=palette[l], edgecolor="white", linewidth=1.5
        )
        for l in labels
    ]
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        fontsize=12,
        frameon=False,
        bbox_to_anchor=(0.5, 0.0),
    )

    fig.suptitle(
        "Wall Load-Bearing Partition by Polygenic Score",
        fontsize=18,
        weight="bold",
        y=1.02,
        color=ink,
    )
    plt.tight_layout(rect=[0, 0.06, 1, 0.95])
    plt.show()
    return fig


def plot_time_step_convergence(dt_list, params):
    """
    Plots the convergence of key variables (systolic stretch, collagen density, and active TGF-β)

    Parameters:
    dt_list (list):
        A list of time step sizes to test for convergence.
    params (dict):
        ArterialParameters / the dictionary containing the parameters for the aneurysm simulation. 
    """
    final_lambda_sys_array = []
    final_collagen_array = []
    final_tgf_array = []

    for dt in dt_list:
        results = model.simulate_aneurysm(params, dt=dt)
        final_lambda_sys_array.append(results["lambda_sys"][-1])
        final_collagen_array.append(results["collagen_ad"][-1])
        final_tgf_array.append(results["active_tgf_beta"][-1])
        print(f"dt={dt:.4f} done")

    # Ratio of current step to previous step
    def step_ratios(arr):
        arr = np.array(arr)
        return arr[1:] / arr[:-1]

    ratio_lam = step_ratios(final_lambda_sys_array)
    ratio_col = step_ratios(final_collagen_array)
    ratio_tgf = step_ratios(final_tgf_array)
    dt_plot = dt_list[1:]  # one shorter than dt_list

    converged_indices = [
        i
        for i in range(len(dt_plot))
        if abs(ratio_lam[i] - 1.0) < 0.01
        and abs(ratio_col[i] - 1.0) < 0.01
        and abs(ratio_tgf[i] - 1.0) < 0.01
    ]
    optimal_dt = dt_plot[converged_indices[0]] if converged_indices else None

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(
        dt_plot,
        ratio_lam,
        marker="o",
        linestyle="-",
        color="black",
        linewidth=2,
        markersize=6,
        label=r"$\lambda_{sys}$",
    )
    ax.plot(
        dt_plot,
        ratio_col,
        marker="x",
        linestyle="-",
        color="#0072B2",
        linewidth=2,
        markersize=6,
        label="Collagen (Adventitia)",
    )
    ax.plot(
        dt_plot,
        ratio_tgf,
        marker=".",
        linestyle="-",
        color="#D55E00",
        linewidth=2,
        markersize=8,
        label="Active TGF-β",
    )

    ax.axhline(
        1.0,
        color="grey",
        linestyle=":",
        linewidth=1.5,
        alpha=0.6,
        label="Ratio = 1 (converged)",
    )

    if optimal_dt is not None:
        ax.axvline(
            optimal_dt,
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Convergence threshold (dt={optimal_dt:.4f} yr)",
        )
        print(f"All variables converged at dt = {optimal_dt:.4f} years")
    else:
        print("No dt achieves convergence within tolerance across all variables.")

    ax.set_xlabel("Time Step (years)", fontsize=20, weight="bold")
    ax.set_ylabel("Ratio: f(dt) / f(dt - 1)", fontsize=18, weight="bold")
    ax.set_title("Time Step Convergence Study", fontsize=22, weight="bold", pad=15)
    ax.legend(fontsize=14, frameon=True, edgecolor="black")
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.show()
    return fig