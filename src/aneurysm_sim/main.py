from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model.model import simulate_arterial_stress_and_pressure, simulate_elastin_degredation, simulate_aneurysm
from aneurysm_sim.model import plots

def main():
    params = ArterialParameters()
    results = simulate_arterial_stress_and_pressure(params)
    # plots.plot_pressure_vs_stretch(results)
    # plots.plot_pressure_vs_diameter(results)

    # Test elastin degradation simulation
    degredation_factors = [1.0, 0.8, 0.6, 0.4, 0.2]
    degredation_results = simulate_elastin_degredation(params, degredation_factors)
    # plots.plot_elastin_degradation(degredation_results)

    # Plot normalised densities
    results_tt = simulate_aneurysm(params, genotype="TT")
    plots.plot_normalised_densities(results_tt)
    results_tc = simulate_aneurysm(params, genotype="TC")
    plots.plot_normalised_densities(results_tc)
    results_cc = simulate_aneurysm(params, genotype="CC")
    plots.plot_normalised_densities(results_cc)
    # plots.plot_systolic_pressure_over_time(results_tt, results_tc, results_cc)

