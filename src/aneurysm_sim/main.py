from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model.model import simulate_arterial_stress_and_pressure, simulate_elastin_degredation, simulate_aneurysm
from aneurysm_sim.model import plots

def main():
    params = ArterialParameters()
    results = simulate_arterial_stress_and_pressure(params)
    # plots.plot_pressure_vs_stretch(results)
    # plots.plot_pressure_vs_diameter(results)

    # Test elastin degradation simulation
    # degredation_factors = [1.0, 0.8, 0.6, 0.4, 0.2]
    # degredation_results = simulate_elastin_degredation(params, degredation_factors)
    # plots.plot_elastin_degradation(degredation_results)

    # # Plot normalised densities
    # results_main = simulate_aneurysm(params)
    # plots.plot_normalised_densities(results_main, title="Normalised Densities Over Time")
    # plots.plot_max_collagen_stretch(results_main)
    results_tt = {
        "No Treatment": simulate_aneurysm(params, genotype="TT", treatment=False),
        "Treatment": simulate_aneurysm(params, genotype="TT", treatment=True),
    }
    results_tc = {
        "No Treatment": simulate_aneurysm(params, genotype="TC", treatment=False),
        "Treatment": simulate_aneurysm(params, genotype="TC", treatment=True),
    }
    results_cc = {
        "No Treatment": simulate_aneurysm(params, genotype="CC", treatment=False),
        "Treatment": simulate_aneurysm(params, genotype="CC", treatment=True),
    }       
    plots.plot_auc_bars_by_genotype(
        results_tt, results_tc, results_cc)
    
    results_tt = simulate_aneurysm(params, genotype="TT")
    # plots.plot_normalised_densities(results_tt)
    results_tc = simulate_aneurysm(params, genotype="TC")
    # plots.plot_normalised_densities(results_tc)
    results_cc = simulate_aneurysm(params, genotype="CC")
    plots.plot_normalised_densities_by_genotype(results_tt, results_tc, results_cc)
    # plots.plot_normalised_densities(results_cc)
    plots.plot_systolic_stretch_over_time(results_tt, results_tc, results_cc)
    plots.plot_tgf_vs_fibroblasts(results_tt, results_tc, results_cc)
    results_treat = simulate_aneurysm(params, treatment=True)
    plots.plot_normalised_densities(results_treat)
    results_treat_tt = simulate_aneurysm(params, genotype="TT", treatment=True)
    results_treat_tc = simulate_aneurysm(params, genotype="TC", treatment=True)
    results_treat_cc = simulate_aneurysm(params, genotype="CC", treatment=True)
    plots.plot_systolic_stretch_over_time(results_treat_tt, results_treat_tc, results_treat_cc)
    plots.plot_normalised_densities_by_genotype(results_treat_tt, results_treat_tc, results_treat_cc)
    plots.plot_tgf_vs_fibroblasts(results_treat_tt, results_treat_tc, results_treat_cc)


    # t_treat_list = [45, 50, 55, 60, 65]

    # results_tt_list = []
    # results_tc_list = []
    # results_cc_list = []

    # for t_treat in t_treat_list:
    #     params = ArterialParameters()
    #     params.t_treat = t_treat
    #     results_tt_list.append(simulate_aneurysm(params, genotype="TT", treatment=True))

    #     params = ArterialParameters()
    #     params.t_treat = t_treat
    #     results_tc_list.append(simulate_aneurysm(params, genotype="TC", treatment=True))

    #     params = ArterialParameters()
    #     params.t_treat = t_treat
    #     results_cc_list.append(simulate_aneurysm(params, genotype="CC", treatment=True))

    # plots.plot_diameter_treatment_times(results_tt_list, results_tc_list, results_cc_list, t_treat_list)