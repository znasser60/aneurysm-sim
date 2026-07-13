# Intracranial Aneurysm Simulation (Finite Element Model)

This repository contains a Python package for simulating an intracranial aneurysm based on varying genotypes and polygenic scores. The project includes both a 1D chemo-mechano-biological model based in Constrained Mixture Theory and a 3D finite element (FEBio) model with a focal immune cell infiltration and visible aneurysm bulge formation. This was created for the purpose of fulfilling the final Master's thesis project for Computational Science program at the University of Amsterdam [July 2026].

## Setup and Installation

The package requires Python 3.8 or higher. 

### 1. Clone the Repository
Clone the repository locally:

**Mac and Windows:**
```bash
git clone [https://github.com/znasser60/aneurysm-sim.git](https://github.com/znasser60/aneurysm-sim.git)
cd aneurysm-sim
```

### 2. Set up the Virtual Environment
Create and activate a virtual environment to manage dependencies:

**Mac:**
```bash
python3 -m venv .venv
source .venv/bin/activate
```

**Windows:**
```bash
python -m venv .venv
.venv\Scripts\activate
```

### 3. Install Dependencies
Install the package and its dependencies:

**Mac and Windows:**
```bash
pip install -e .
```

---

## 1D Model Simulation

The 1D model simulates individual cell and protein density evolution and the resulting systolic stretch during an intracranial aneurysm over 50 years for a hypothetical patient. You can run the model using the installed command-line tool `aneurysm-sim` or by calling the module directly with `python -m aneurysm_sim`. 

The CLI works in two primary modes: `general` (population-level comparisons) and `patient` (individual patient simulations).

Two pieces of genetic information can be varied: 
1. TGF-B genotype: A genotype that dictates the baseline production levels of TGF-B (options: TT, TC, CC)
2. vSMC polygenic score: A score that indicates the number of risk alleles that impact vSMC volume fractions (options: 0-4; a higher value generally suggests less smooth muscle cells in the artery) 

### General Mode
Runs population-level comparisons across the three TGF-beta genotypes (TT, TC, CC) and/or five polygenic risk scores (0-4): 

**CLI Command:**
```bash
aneurysm-sim --mode general --plot <plot_name>
```

**Available general plots**:
*   `density_by_genotype`: Normalised constituent mass density evolution across genotypes, no treatment.
*   `density_by_geno_treat`: Normalised constituent mass density evolution across genotypes, with treatment.
*   `density_by_score`: Normalised constituent density evolution across polygenic scores, no treatment. 
*   `density_by_score_treat`: Normalised constituent density evolution across polygenic scores, with treatment.
*   `stretch_by_genotype`: Systolic stretch trajectories over time across genotypes.
*   `stretch_by_geno_treat`: Side-by-side systolic stretch trajectories by genotype, without vs. with treatment.
*   `stretch_by_score`: Systolic stretch trajectories by polygenic score (plots sensitivity analysis with 100 runs).
*   `stretch_by_score_treat`: Systolic stretch trajectories by polygenic score, shown side-by-side without vs. with treatment (plots sensitivity analysis with 100 runs).
*   `load_bearing`: Collagen vs. smooth-muscle load-fraction partitioning over time.
*   `time_convergence`: Numerical convergence of the model output as the integration time step is refined.

### Patient Mode
Runs simulations for specific individual patients using genotype and polygenic score values looked up from the provided dataset by patient ID. This mode requires the `--id` flag. 

See `src/aneurysm_sim/data/syn_patient_data.csv` to choose a hypothetical patient with specified genetic information. 

**CLI Command:**
```bash
aneurysm-sim --mode patient --id <patient_id> --plot <plot_name>
```
*(NOTE: You can pass multiple IDs and multiple plots by separating them with spaces.)*

**Available patient plots**:
*   `pressure_vs_stretch`: Modelled arterial pressure-stretch relationship.
*   `pressure_vs_diameter`: Modelled arterial pressure-diameter relationship.
*   `stretch_vs_stress`: Modelled stretch-stress relationship.
*   `att_dist`: Evolving distribution of collagen fibre attachment stretch.
*   `rec_dist`: Evolving distribution of collagen fibre recruitment stretch.
*   `density`: Normalised constituent density evolution, no treatment.
*   `density_treat`: Normalised constituent density evolution, with treatment.
*   `stretch`: Systolic stretch over time.
*   `stretch_treat`: Systolic stretch over time, with treatment.

---

## 3D FEBio Model Simulation

The 3D model workflow uses FEBio, the open source nonlinear finite element solver designed for computational biomechanics. This pipeline generates patient-specific FEBio `.feb` input files based on user input of genotype and polygenic score. Fitted Holzapfel-Gasser-Ogden (HGO) material parameters are placed in these files and launched in the FEBio software to run a custom aneurysm plugin `.dylib`. 

**Prerequisite:** Running the simulations requires FEBioStudio to be installed on your machine. The runner script looks for the Mac FEBio executable at `/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio4` by default.

### 1. Run the FEBio Simulations
You can execute the FEBio runner script to run a batch of simulations in parallel or with a queue. The runner accepts `--genotypes` and `--scores` flags to define which simulations will be ran. 

**Mac and Windows:**
```bash
python src/aneurysm_sim/febio/febio_runner.py --genotypes TT TC CC --scores 0 4
```

### 2. Plot the FEBio Results
Once the simulations are complete, the resulting `.xplt` data is extracted to CSV files containing effective stress and x-displacement data. You can visualize these results using the FEBio plotter module.

**Mac and Windows:**
```bash
python src/aneurysm_sim/febio/febio_plotter.py
```

This will generate the following plots:
*   Relative x displacement over time for each genotype.

*   Relative effective stress change over time for each genotype.
*   Relative x displacement over time for each polygenic score.
*   Relative effective stress change over time for each polygenic score.

Final 3D aneurysm simulations and data files can be found in `src/aneurysm_sim/febio/`. 
