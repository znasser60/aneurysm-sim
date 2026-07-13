import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import model, plots

DATA_PATH = "src/aneurysm_sim/data/syn_patient_data.csv"


def run_general_mode(plot_names):
    GENOTYPES = ["TT", "TC", "CC"]
    SCORES = range(5)
    REDUCED_SCORES = [0, 2, 4]
    SMC_RANGE = np.linspace(0.45, 0.85, 30)
    TGF_RANGE = np.linspace(0.65, 1.25, 30)
    Z = np.zeros((len(TGF_RANGE), len(SMC_RANGE)))

    for plot_name in plot_names:
        if plot_name == "density_by_genotype":
            results = {
                g: model.simulate_aneurysm(ArterialParameters(genotype=g))
                for g in GENOTYPES
            }
            plots.plot_normalised_densities_by_genotype(
                results["TT"], results["TC"], results["CC"]
            )

        elif plot_name == "density_by_geno_treat":
            results = {
                g: model.simulate_aneurysm(
                    ArterialParameters(genotype=g), treatment=True
                )
                for g in GENOTYPES
            }
            plots.plot_normalised_densities_by_genotype(
                results["TT"], results["TC"], results["CC"]
            )

        elif plot_name == "density_by_score":
            results = {
                s: model.simulate_aneurysm(ArterialParameters(polygenic_score=s))
                for s in REDUCED_SCORES
            }
            plots.plot_normalised_densities_by_score(
                *[results[s] for s in REDUCED_SCORES]
            )

        elif plot_name == "density_by_score_treat":
            results = {
                s: model.simulate_aneurysm(
                    ArterialParameters(polygenic_score=s), treatment=True
                )
                for s in REDUCED_SCORES
            }
            plots.plot_normalised_densities_by_score(
                *[results[s] for s in REDUCED_SCORES]
            )

        elif plot_name == "stretch_by_genotype":
            results = {
                g: model.simulate_aneurysm(ArterialParameters(genotype=g))
                for g in GENOTYPES
            }
            plots.plot_stretch_by_genotype(results["TT"], results["TC"], results["CC"])

        elif plot_name == "stretch_by_geno_treat":
            results = {
                g: model.simulate_aneurysm(ArterialParameters(genotype=g))
                for g in GENOTYPES
            }
            results_treat = {
                g: model.simulate_aneurysm(
                    ArterialParameters(genotype=g), treatment=True
                )
                for g in GENOTYPES
            }
            plots.plot_stretch_by_geno_treat(
                results["TT"],
                results["TC"],
                results["CC"],
                results_treat["TT"],
                results_treat["TC"],
                results_treat["CC"],
            )

        elif plot_name == "stretch_by_score":
            results_batch = {}
            for s in SCORES:
                p = ArterialParameters(polygenic_score=s)
                results_batch[s] = model.simulate_aneurysm_batch_smc(
                    p, treatment=False, n_vals=100
                )
            plots.plot_stretch_by_score(results_batch)

        elif plot_name == "stretch_by_score_treat":
            results_notreat, results_treat = {}, {}
            for s in SCORES:
                p = ArterialParameters(polygenic_score=s)
                results_notreat[s] = model.simulate_aneurysm_batch_smc(
                    p, treatment=False, n_vals=1
                )
                results_treat[s] = model.simulate_aneurysm_batch_smc(
                    p, treatment=True, n_vals=1
                )
            plots.plot_stretch_by_score(results_notreat, results_treat)

        elif plot_name == "load_bearing":
            params_dict = {s: ArterialParameters(polygenic_score=s) for s in SCORES}
            plots.plot_load_bearing(params_dict)

        elif plot_name == "time_convergence":
            dt_list = np.linspace(0.05, 0.0005, 50)
            plots.plot_time_step_convergence(dt_list, ArterialParameters())

        else:
            print(
                f"Unknown plot '{plot_name}' for general mode. Available: density_by_genotype, density_by_geno_treat, density_by_score"
            )


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
        patient_stretch_results = model.simulate_arterial_stress_and_pressure(
            patient_params
        )

        plot_registry = {
            "pressure_vs_stretch": lambda: plots.plot_pressure_vs_stretch(
                patient_stretch_results
            ),
            "pressure_vs_diameter": lambda: plots.plot_pressure_vs_diameter(
                patient_stretch_results
            ),
            "stretch_vs_stress": lambda: plots.plot_stretch_vs_stress(
                patient_stretch_results
            ),
            "att_dist": lambda: plots.plot_stretch_evolution_dist(
                patient_results, stretch_type="attachment"
            ),
            "rec_dist": lambda: plots.plot_stretch_evolution_dist(
                patient_results, stretch_type="recruitment"
            ),
            "density": lambda: [
                plots.plot_normalised_densities(patient_results, legend=True),
                plt.show(),
            ],
            "density_treat": lambda: [
                plots.plot_normalised_densities(patient_results_treat, legend=True),
                plt.show(),
            ],
            "stretch": lambda: [
                plots.plot_stretch(patient_results),
                plt.show(),
            ],
            "stretch_treat": lambda: [
                plots.plot_stretch(patient_results_treat),
                plt.show(),
            ],
        }

        for plot_name in plot_names:
            if plot_name in plot_registry:
                plot_registry[plot_name]()
            else:
                print(
                    f"Unknown plot '{plot_name}' for patient mode. Available: {', '.join(plot_registry)}"
                )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Genotype- and polygenic-score-dependent chemo-bio-mechanical"
            "simulation of intracranial aneurysm growth. Generates comparison plots of "
            "constituent density evolution, systolic stretch trajectories, load-bearing "
            "epochs, and numerical convergence, either aggregated across all "
            "genotype/score cohorts ('general' mode) or for individual patients from a "
            "supplied dataset ('patient' mode)."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["patient", "general"],
        default="general",
        help=(
            "Simulation mode. 'general': runs population-level comparisons across the three "
            "TGF-beta genotypes (TT/TC/CC) and/or the five polygenic risk scores (0-4), each "
            "simulated from ArterialParameters defaults. 'patient': runs simulations for one "
            "or more individual patients, using genotype and polygenic score values looked up "
            "from the CSV dataset by patient ID."
        ),
    )
    parser.add_argument(
        "--id",
        type=int,
        nargs="+",
        help=(
            "One or more patient IDs to simulate, matched against the 'ID' column of the "
            "patient dataset. Required when --mode patient is selected, ignored otherwise."
        ),
    )
    parser.add_argument(
        "--plot",
        type=str,
        nargs="+",
        help=(
            "One or more plot names to generate.\n"
            "\n"
            "General mode (population-level, by genotype or polygenic score):\n"
            "  density_by_genotype    Normalised constituent mass density evolution, "
            "compared across TT/TC/CC genotypes, no treatment.\n"
            "  density_by_geno_treat  As above, with treatment intervention applied.\n"
            "  density_by_score       Normalised constituent density evolution, compared "
            "across polygenic scores 0-4, no treatment.\n"
            "  density_by_score_treat As above, with treatment intervention applied.\n"
            "  stretch_by_genotype    Systolic stretch (lambda_sys) trajectories over time, "
            "compared across TT/TC/CC genotypes.\n"
            "  stretch_by_geno_treat  Side-by-side systolic stretch trajectories by genotype, "
            "without vs. with treatment intervention.\n"
            "  stretch_by_score       Systolic stretch trajectories by polygenic score, with "
            "shaded min/max envelopes from a Gaussian vSMC-fraction sensitivity sweep "
            "(deterministic, not Monte Carlo) about each score's mean.\n"
            "  stretch_by_score_treat As above, shown side-by-side without vs. with "
            "treatment intervention.\n"
            "  load_bearing_epochs    Collagen vs. smooth-muscle load-fraction partitioning "
            "over time, contrasting the lowest (0) and highest (4) polygenic score.\n"
            "  time_convergence       Numerical convergence of the model output as the "
            "integration time step dt is refined, for validation of temporal discretisation.\n"
            "\n"
            "Patient mode (individual patient, by ID):\n"
            "  pressure_vs_stretch    Modelled arterial pressure-stretch relationship for "
            "the patient's fitted material parameters.\n"
            "  pressure_vs_diameter   Modelled arterial pressure-diameter relationship.\n"
            "  stretch_vs_stress      Modelled stretch-stress (constitutive) relationship.\n"
            "  att_dist               Evolving distribution of collagen fibre attachment "
            "stretch over the simulation.\n"
            "  rec_dist               Evolving distribution of collagen fibre recruitment "
            "stretch over the simulation.\n"
            "  density                Normalised constituent density evolution for the "
            "patient, no treatment.\n"
            "  density_treat          As above, with treatment intervention applied.",
        ),
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
