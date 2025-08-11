import argparse

from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import model
from aneurysm_sim.model import plots

def main():
    parser = argparse.ArgumentParser(description="Run aneurysm simulations and plots.")
    parser.add_argument(
        "--plot",
        type=str,
        nargs="+",
        help="Which plot(s) to run. Options: stretch_vs_stress, att_dist, rec_dist, "
             "norm_density, stretch_time, norm_density_treat, stretch_time_treat, auc_bars"
    )
    args = parser.parse_args()

    params = ArterialParameters()

    results = model.simulate_arterial_stress_and_pressure(params)
    results_tt = model.simulate_aneurysm(params, genotype="TT")
    results_tc = model.simulate_aneurysm(params, genotype="TC")
    results_cc = model.simulate_aneurysm(params, genotype="CC")
    results_treat_tt = model.simulate_aneurysm(params, genotype="TT", treatment=True)
    results_treat_tc = model.simulate_aneurysm(params, genotype="TC", treatment=True)
    results_treat_cc = model.simulate_aneurysm(params, genotype="CC", treatment=True)

    results_tt_dict = {
        "No Treatment": results_tt,
        "Treatment": results_treat_tt,
    }
    results_tc_dict = {
        "No Treatment": results_tc,
        "Treatment": results_treat_tc,
    }
    results_cc_dict = {
        "No Treatment": results_cc,
        "Treatment": results_treat_cc,
    }

    # Run only the selected plots
    if not args.plot:
        print("No plot selected. Use --plot to specify one or more plots.")
        return

    for plot_name in args.plot:
        if plot_name == "stretch_vs_stress":
            plots.plot_stretch_vs_stress(results)
        elif plot_name == "att_dist":
            plots.plot_att_dist(params.c_att_min_ad, params.c_att_mod_ad, params.c_att_max_ad)
        elif plot_name == "rec_dist":
            plots.plot_rec_dist(params.c_rec_min_ad, params.c_rec_mod_ad, params.c_rec_max_ad)
        elif plot_name == "norm_density":
            plots.plot_normalised_densities_by_genotype(results_tt, results_tc, results_cc)
        elif plot_name == "sys_stretch":
            plots.plot_systolic_stretch_over_time(results_tt, results_tc, results_cc)
        elif plot_name == "norm_density_treat":
            plots.plot_normalised_densities_by_genotype(results_treat_tt, results_treat_tc, results_treat_cc)
        elif plot_name == "sys_stretch_treat":
            plots.plot_systolic_stretch_over_time(results_treat_tt, results_treat_tc, results_treat_cc)
        elif plot_name == "auc_bars":
            plots.plot_auc_bars_by_genotype(results_tt_dict, results_tc_dict, results_cc_dict)
        else:
            print(f"Unknown plot: {plot_name}")