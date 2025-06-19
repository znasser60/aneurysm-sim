from aneurysm_sim.model.model import simulate_arterial_stress_and_pressure
from aneurysm_sim.model import plots

def main():
    results = simulate_arterial_stress_and_pressure()
    plots.plot_pressure_vs_stretch(results)
    plots.plot_pressure_vs_diameter(results)

if __name__ == "__main__":
    main()