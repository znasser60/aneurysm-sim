from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model.model import simulate_arterial_stress_and_pressure, simulate_aneurysm
from aneurysm_sim.model import plots

def main():
    params = ArterialParameters()
    results = simulate_arterial_stress_and_pressure(params)
    plots.plot_stretch_vs_stress(results)
    plots.plot_att_dist(params.c_att_min_ad, params.c_att_mod_ad, params.c_att_max_ad)
    plots.plot_rec_dist(params.c_rec_min_ad, params.c_rec_mod_ad, params.c_rec_max_ad)
    results_tt = simulate_aneurysm(params, genotype="TT")
    results_tc = simulate_aneurysm(params, genotype="TC")
    results_cc = simulate_aneurysm(params, genotype="CC")
    plots.plot_normalised_densities_by_genotype(results_tt, results_tc, results_cc)
    plots.plot_systolic_stretch_over_time(results_tt, results_tc, results_cc)
    results_treat_tt = simulate_aneurysm(params, genotype="TT", treatment=True)
    results_treat_tc = simulate_aneurysm(params, genotype="TC", treatment=True)
    results_treat_cc = simulate_aneurysm(params, genotype="CC", treatment=True)
    plots.plot_normalised_densities_by_genotype(results_treat_tt, results_treat_tc, results_treat_cc)
    plots.plot_systolic_stretch_over_time(results_treat_tt, results_treat_tc, results_treat_cc)

    results_tt_dict = {
        "No Treatment": simulate_aneurysm(params, genotype="TT", treatment=False),
        "Treatment": simulate_aneurysm(params, genotype="TT", treatment=True),
    }
    results_tc_dict = {
        "No Treatment": simulate_aneurysm(params, genotype="TC", treatment=False),
        "Treatment": simulate_aneurysm(params, genotype="TC", treatment=True),
    }
    results_cc_dict = {
        "No Treatment": simulate_aneurysm(params, genotype="CC", treatment=False),
        "Treatment": simulate_aneurysm(params, genotype="CC", treatment=True),
    }       
    plots.plot_auc_bars_by_genotype(
        results_tt_dict, results_tc_dict, results_cc_dict)
