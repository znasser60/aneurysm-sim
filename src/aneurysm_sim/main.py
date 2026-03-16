import pandas as pd
import matplotlib.pyplot as plt
import argparse

from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import model, plots

DATA_PATH = "src/aneurysm_sim/data/syn_patient_data.csv"


def run_general_mode(plot_names):
    """Run only the requested general comparison simulations and plots."""
    GENOTYPES = ["TT", "TC", "CC"]
    SCORES = range(5)

    for plot_name in plot_names:
        if plot_name == "density_by_genotype":
            results = {g: model.simulate_aneurysm(ArterialParameters(genotype=g)) for g in GENOTYPES}
            plots.plot_normalised_densities_by_genotype(results["TT"], results["TC"], results["CC"])

        elif plot_name == "density_by_geno_treat":
            results = {g: model.simulate_aneurysm(ArterialParameters(genotype=g), treatment=True) for g in GENOTYPES}
            plots.plot_normalised_densities_by_genotype(results["TT"], results["TC"], results["CC"])

        elif plot_name == "density_by_score":
            results = {s: model.simulate_aneurysm(ArterialParameters(polygenic_score=s)) for s in SCORES}
            plots.plot_normalised_densities_by_score(*[results[s] for s in SCORES])
        
        elif plot_name == "density_by_score_treat":
            results = {s: model.simulate_aneurysm(ArterialParameters(polygenic_score=s), treatment=True) for s in SCORES}
            plots.plot_normalised_densities_by_score(*[results[s] for s in SCORES])

        elif plot_name == "stretch_by_genotype":
            results = {g: model.simulate_aneurysm(ArterialParameters(genotype=g)) for g in GENOTYPES}
            plots.plot_stretch_by_genotype(results["TT"], results["TC"], results["CC"])
        
        elif plot_name == "stretch_by_geno_treat":
            results = {g: model.simulate_aneurysm(ArterialParameters(genotype=g), treatment=True) for g in GENOTYPES}
            plots.plot_stretch_by_genotype(results["TT"], results["TC"], results["CC"])

        elif plot_name == "stretch_by_score":
            results = {s: model.simulate_aneurysm(ArterialParameters(polygenic_score=s)) for s in SCORES}
            plots.plot_stretch_by_score(*[results[s] for s in SCORES])

        elif plot_name == "stretch_by_score_treat":
            results = {s: model.simulate_aneurysm(ArterialParameters(polygenic_score=s), treatment=True) for s in SCORES}
            plots.plot_stretch_by_score(*[results[s] for s in SCORES])

        else:
            print(f"Unknown plot '{plot_name}' for general mode. Available: density_by_genotype, density_by_geno_treat, density_by_score")


def run_patient_mode(patient_data, patient_ids, plot_names):
    """Run only patient-specific plots — no population comparisons."""
    for patient_id in patient_ids:
        patient_row = patient_data[patient_data["ID"] == patient_id]
        if patient_row.empty:
            print(f"Patient ID {patient_id} not found in the dataset.")
            continue

        patient_params = ArterialParameters(
            genotype=patient_row["TGF Genotype"].values[0],
            polygenic_score=patient_row["Polygenic Score"].values[0],
        )
        
        patient_results = model.simulate_aneurysm(patient_params)
        patient_results_treat = model.simulate_aneurysm(patient_params, treatment=True)
        patient_stretch_results = model.simulate_arterial_stress_and_pressure(patient_params)

        plot_registry = {
            "pressure_vs_stretch": lambda: plots.plot_pressure_vs_stretch(patient_stretch_results),
            "pressure_vs_diameter": lambda: plots.plot_pressure_vs_diameter(patient_stretch_results),
            "stretch_vs_stress": lambda: plots.plot_stretch_vs_stress(patient_stretch_results),
            "att_dist":          lambda: plots.plot_att_dist(patient_params.c_att_min_ad, patient_params.c_att_mod_ad, patient_params.c_att_max_ad),
            "rec_dist":          lambda: plots.plot_rec_dist(patient_params.c_rec_min_ad, patient_params.c_rec_mod_ad, patient_params.c_rec_max_ad),
            "density":           lambda: [plots.plot_normalised_densities(patient_results, legend=True), plt.show()],
            "density_treat":     lambda: [plots.plot_normalised_densities(patient_results_treat, legend=True), plt.show()],
        }

        for plot_name in plot_names:
            if plot_name in plot_registry:
                plot_registry[plot_name]()
            else:
                print(f"Unknown plot '{plot_name}' for patient mode. Available: {', '.join(plot_registry)}")


def main():
    parser = argparse.ArgumentParser(description="Run aneurysm simulations and plots.")
    parser.add_argument(
        "--mode",
        type=str,
        choices=["patient", "general"],
        default="general",
        help="'patient' for individual patient data only, 'general' for population-level comparisons."
    )
    parser.add_argument(
        "--id",
        type=int,
        nargs="+",
        help="Patient ID(s) from the CSV. Required for patient mode."
    )
    parser.add_argument(
        "--plot",
        type=str,
        nargs="+",
        help=(
            "Plot name(s) to generate.\n"
            "  General mode: density_by_genotype, density_by_geno_treat, density_by_score\n"
            "  Patient mode: stretch_vs_stress, att_dist, rec_dist, density, density_treat"
        )
    )
    args = parser.parse_args()

    if not args.plot:
        print("No plot selected. Use --plot to specify one or more plots.")
        return

    if args.mode == "patient" and not args.id:
        print("Patient mode requires --id. Please specify one or more patient IDs.")
        return

    if args.mode == "general":
        run_general_mode(args.plot)
    else:
        patient_data = pd.read_csv(DATA_PATH)
        run_patient_mode(patient_data, args.id, args.plot)
