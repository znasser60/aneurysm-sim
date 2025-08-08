from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import model
from aneurysm_sim.model import plots

def main():
    params = ArterialParameters()
    results = model.simulate_arterial_stress_and_pressure(params)
    plots.plot_stretch_vs_stress(results)
    plots.plot_att_dist(params.c_att_min_ad, params.c_att_mod_ad, params.c_att_max_ad)
    plots.plot_rec_dist(params.c_rec_min_ad, params.c_rec_mod_ad, params.c_rec_max_ad)
    results_tt = model.simulate_aneurysm(params, genotype="TT")
    results_tc = model.simulate_aneurysm(params, genotype="TC")
    results_cc = model.simulate_aneurysm(params, genotype="CC")
    print("Final Lambda Values: " \
    "TT: {}, TC: {}, CC: {}".format(
        results_tt["final_lambda_sys"], results_tc["final_lambda_sys"], results_cc["final_lambda_sys"]
    ))
    plots.plot_normalised_densities_by_genotype(results_tt, results_tc, results_cc)
    plots.plot_systolic_stretch_over_time(results_tt, results_tc, results_cc)
    results_treat_tt = model.simulate_aneurysm(params, genotype="TT", treatment=True)
    results_treat_tc = model.simulate_aneurysm(params, genotype="TC", treatment=True)
    results_treat_cc = model.simulate_aneurysm(params, genotype="CC", treatment=True)
    print("Final Lambda Values With Treatment: " \
    "TT: {}, TC: {}, CC: {}".format(
        results_treat_tt["final_lambda_sys"], results_treat_tc["final_lambda_sys"], results_treat_cc["final_lambda_sys"]
    ))
    plots.plot_normalised_densities_by_genotype(results_treat_tt, results_treat_tc, results_treat_cc)
    plots.plot_systolic_stretch_over_time(results_treat_tt, results_treat_tc, results_treat_cc)

    results_tt_dict = {
        "No Treatment": model.simulate_aneurysm(params, genotype="TT", treatment=False),
        "Treatment": model.simulate_aneurysm(params, genotype="TT", treatment=True),
    }
    results_tc_dict = {
        "No Treatment": model.simulate_aneurysm(params, genotype="TC", treatment=False),
        "Treatment": model.simulate_aneurysm(params, genotype="TC", treatment=True),
    }
    results_cc_dict = {
        "No Treatment": model.simulate_aneurysm(params, genotype="CC", treatment=False),
        "Treatment": model.simulate_aneurysm(params, genotype="CC", treatment=True),
    }
    plots.plot_auc_bars_by_genotype(
        results_tt_dict, results_tc_dict, results_cc_dict)
