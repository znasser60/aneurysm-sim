"""This module generates patient-specific FEBio .feb inputs and run the aneurysm plugin.
 
Injects genotype, polygenic score, and fitted HGO material parameters into a
template .feb file, then launches FEBio with the custom plugin over a batch
of genotype/score setups.
"""
import xml.etree.ElementTree as ET
import subprocess
import os
import time
import argparse
import sys

from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model.febio_solver import (
    fit_media_stress_hgo_params,
    fit_adventitia_stress_hgo_params,
)

current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, "../../"))
sys.path.insert(0, src_dir)


class FebioRunner:
    """
    Class to generate patient-specific FEBio .feb files and run simulations with the aneurysm plugin.

    Attributes:
    template_path (str):
        Path to the template .feb file.
    results_dir (str):
        Directory to save the generated .feb files and simulation results.
    febio_executable (str):
        Path to the FEBio executable.
    plugin_path (str):
        Path to the compiled aneurysm plugin shared library.
    
    Methods:
    _generate_patient_feb(params, output_path):
        Generates a patient-specific .feb file with the given parameters.
    run_simulations(batch_setups, max_concurrent=3, cores_per_sim=2):
        Runs a batch of simulations in parallel or in a queue with the specified setups.
    """
    def __init__(
        self,
        template_path="cylinder_with_plugin.feb",
        results_dir="src/aneurysm_sim/febio/results",
    ):
        self.template_path = template_path
        if results_dir is None:
            results_dir = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "results"
            )
        self.results_dir = results_dir
        os.makedirs(self.results_dir, exist_ok=True)

        self.febio_executable = (
            "/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio4"
        )
        current_dir = os.path.dirname(os.path.abspath(__file__))
        self.plugin_path = os.path.join(current_dir, "build", "AneurysmPlugin.dylib")

    def _generate_patient_feb(self, params: ArterialParameters, output_path: str):
        """
        Loads the template XML, modifies the genotype, score, and HGO parameters,
        and saves it to a new file.
        """
        tree = ET.parse(self.template_path)
        root = tree.getroot()

        geno_map = {"TT": 0, "TC": 1, "CC": 2}
        geno_int = geno_map.get(params.genotype.upper(), 0)
        score_int = (
            int(params.polygenic_score) if params.polygenic_score is not None else 0
        )
        for mat_name in ["Media", "Adventitia"]:
            mat = root.find(f".//material[@name='{mat_name}']")
            if mat is not None:
                for tag, val in [
                    ("genotype", geno_int),
                    ("polygenic_score", score_int),
                ]:
                    elem = mat.find(tag)
                    if elem is None:
                        elem = ET.SubElement(mat, tag)
                    elem.text = str(val)

        # Get new material parameters based on SMC fraction
        media_c, media_k1, media_k2 = fit_media_stress_hgo_params(params)
        adv_c, adv_k1, adv_k2 = fit_adventitia_stress_hgo_params(params)

        # Update media
        media_mat = root.find(".//material[@name='Media']")
        if media_mat is not None:
            media_mat.find("genotype").text = str(geno_int)
            media_mat.find("polygenic_score").text = str(score_int)

            hgo_media = media_mat.find(".//passive_solid")
            if hgo_media is not None:
                hgo_media.find("c").text = str(media_c)
                hgo_media.find("k1").text = str(media_k1)
                hgo_media.find("k2").text = str(media_k2)

            active_muscle = media_mat.find(".//active_solid")
            if active_muscle is not None:
                smax_pa = params.k_muscle_a * params.vasodil_conc
                smax_mpa = smax_pa / 1e6
                active_muscle.find("smax").text = f"{smax_mpa:.10f}"

        # Update adventitia
        adv_mat = root.find(".//material[@name='Adventitia']")
        if adv_mat is not None:
            hgo_adv = adv_mat.find(".//passive_solid")
            if hgo_adv is not None:
                hgo_adv.find("c").text = str(adv_c)
                hgo_adv.find("k1").text = str(adv_k1)
                hgo_adv.find("k2").text = str(adv_k2)
        results_abs = os.path.abspath(self.results_dir)
        stem = os.path.splitext(os.path.basename(output_path))[0]

        # Output file paths for FEBio to write results
        output_elem = root.find(".//Output")
        if output_elem is not None:
            plotfile = output_elem.find("plotfile")
            if plotfile is not None:
                plotfile.set("file", os.path.join(results_abs, f"{stem}.xplt"))
            logfile = output_elem.find("logfile")
            if logfile is not None:
                logfile.set("file", os.path.join(results_abs, f"{stem}.log"))

        tree.write(output_path, encoding="ISO-8859-1", xml_declaration=True)

    def run_simulations(self, batch_setups: list, max_concurrent=3, cores_per_sim=2):
        """
        Runs a batch of simulations in parallel or in a queue with the specified setups.

        Parameters:
        batch_setups (list):
            List of ArterialParameters instances for each simulation setup.
        max_concurrent (int, optional):
            Maximum number of concurrent simulations to run. Default is 3.
        cores_per_sim (int, optional):
            Number of CPU cores to allocate per simulation. Default is 2. 
            
        NOTE: Do not exceed the total available cores on your machine, it will cause issues.
        """
        print(f"Preparing {len(batch_setups)} simulations...")
        commands = []
        for i, params in enumerate(batch_setups):
            filename = f"sim_{params.genotype}_score{params.polygenic_score}.feb"
            out = os.path.join(self.results_dir, filename)
            self._generate_patient_feb(params, out)
            commands.append(
                (
                    filename,
                    [self.febio_executable, "-i", out, "-import", self.plugin_path],
                )
            )

        active = []
        start_time = time.time()
        print(f"Executing parallel queue (Max {max_concurrent})...")

        for filename, cmd in commands:
            while len(active) >= max_concurrent:
                active = [p for p in active if p.poll() is None]
                time.sleep(0.5)

            print(f"-> Launching: {filename}")
            active.append(
                subprocess.Popen(
                    cmd,
                    stdout=subprocess.DEVNULL,
                )
            )

        for p in active:
            p.wait()

        print(f"\nFinished in {(time.time() - start_time) / 60:.2f} minutes.")


def main():
    parser = argparse.ArgumentParser(description="Run FEBio Aneurysm Simulations.")
    parser.add_argument(
        "--genotypes",
        nargs="+",
        type=str,
        default=["TT", "TC", "CC"],
        help="Which genotypes to run. Example: --genotypes TT CC",
    )
    parser.add_argument(
        "--scores",
        nargs="+",
        type=int,
        default=[0],
        help="Which polygenic scores to run (0-4). Example: --scores 0 4",
    )
    args = parser.parse_args()
    template_file = os.path.join(current_dir, "cylinder_with_plugin.feb")
    runner = FebioRunner(template_path=template_file)
    setups_to_run = []
    for genotype in args.genotypes:
        for score in args.scores:
            params = ArterialParameters(genotype=genotype, polygenic_score=score)
            setups_to_run.append(params)

    if setups_to_run:
        runner.run_simulations(
            batch_setups=setups_to_run, max_concurrent=1, cores_per_sim=2
        )
    else:
        print("No setups defined. Exiting.")


if __name__ == "__main__":
    main()
