# Aneurysm Simulation

Created for the Individual Project course (August 2025). 

This repository contains a Python package for simulating an aneurysm based on varying genotype that determines baseline TGF-β levels. 

## Setup

### 1. Clone the Repository
```bash
git clone https://github.com/znasser60/aneurysm-sim.git
cd aneurysm-sim
```

### 2. Set up the environment
```bash
python3 -m venv .venv
source .venv/bin/activate     # Mac/Linux
# .venv\Scripts\activate      # Windows
```

### 3. Install dependencies 
```bash
pip install -e .
```

### 4. Run the plots from the project root with the ```--plot``` flag.
Available plots: 
- ```stretch_vs_stress```: Plots stress-stretch relationship for medial elastin and collagen, and adventitial collagen
- ```att_dist```: Plots the triangular attachment stretch distribution of the adventitial collgen fibres
- ```rec_dist```: Plots the triangular recruitment stretch distribution of the adventitial collgen fibres
- ```norm_density```: Plots normalised numerical densities of components in the arterial wall over time of simulation _without_ treatment for each genotype
- ```sys_stretch```: Plots systolic stretch over time of simulation _without_ treatment for each genotype
- ```norm_density_treat```: Plots normalised numerical densities of components in the arterial wall over time of simulation _with_ treatment for each genotype
- ```sys_stretch_treat```: Plots systolic stretch over time of simulation _with_ treatment for each genotype
- ```auc_bars```: Plots area under the curve of fibroblasts, adventitial collagen, latent TGF-β, and active TGF-β over the simulation for each genotype 

#### Example: Run a single plot 
```bash
python -m aneurysm_sim --plot norm_density
```

#### Example: Run multiple plots
```bash
python -m aneurysm_sim --plot stretch_vs_stress rec_dist auc_bars
```





