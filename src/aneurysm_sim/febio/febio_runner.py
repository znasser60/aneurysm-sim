import xml.etree.ElementTree as ET
import subprocess
import os

from aneurysm_sim.config.parameters import ArterialParameters
from model.febio_solver import fit_media_stress_hgo_params, fit_adventitia_stress_hgo_params

class FebioRunner:
    def __init__(self, template_path = "cylinder_with_plugin.feb", run_path="patient_sim.feb"):
        self.template_path = template_path
        self.run_path = run_path
        self.results_dir = "results"
        os.makedirs(self.results_dir, exist_ok=True)

    def run_patient_simulation(self, params: ArterialParameters):
        pass