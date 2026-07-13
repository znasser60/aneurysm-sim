"""
This module is a class for plotting FEBio simulation results for aneurysm simulations.
It reads CSV files containing effective stress and x displacement data for different genotypes 
and polygenic scores (see src/data/), directly taken from .xplt data from FEBio simulations. 
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "results")
DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data")


class FEBioPlotter:
    """
    Class for plotting FEBio simulation results for aneurysm simulations.

    Attributes:
    genotypes (list):
        List of genotypes to plot.
    scores (list): 
        List of polygenic scores to plot.
    measurements (list): 
        List of measurements to plot.
    data (dict): 
        Dictionary to store the loaded data for each genotype and score.
    data_path (str): 
        Path to the directory containing the CSV data files.
    
    Methods:
    plot_displacement_by_geno(): 
        Plots relative x displacement over time for each genotype.
    plot_stress_by_geno(): 
        Plots relative effective stress change over time for each genotype
    plot_displacement_by_score(): 
        Plots relative x displacement over time for each polygenic score.
    plot_stress_by_score(): 
        Plots relative effective stress change over time for each polygenic score. 
    """
    def __init__(self, data_path=DATA_PATH):
        self.genotypes = ["TT", "TC", "CC"]
        self.scores = [0, 4]
        self.measurements = ["eff_stress", "x_disp"]
        self.data = {}
        self.data_path = data_path

        for geno in self.genotypes:
            for score in self.scores:
                dfs = []
                for meas in self.measurements:
                    file_path = os.path.join(data_path, f"{meas}_{geno}_{score}.csv")
                    if not os.path.exists(file_path):
                        continue
                    df = pd.read_csv(file_path, sep=r"\s+").rename(
                        columns={"E1095": meas}
                    )
                    dfs.append(df)
                if dfs:
                    merged = dfs[0]
                    for df in dfs[1:]:
                        merged = merged.merge(df, on="x", how="outer")
                    self.data[(geno, score)] = merged.sort_values("x").reset_index(
                        drop=True
                    )

    def plot_displacement_by_geno(self):
        """Plot relative x displacement over time for each genotype."""
        plt.figure()
        for geno in self.genotypes:
            df = self.data.get((geno, 0))
            if df is not None:
                time_years = (df["x"] - 1.1) * 90
                idx = (df["x"] - 1.1).abs().argmin()
                homeostasis_val = df["x_disp"].iloc[idx]
                disp = (df["x_disp"] - homeostasis_val) / homeostasis_val
                plt.plot(time_years, disp, label=f"{geno} score 0")
        plt.xlabel("Time (years)", fontsize=20)
        plt.ylabel("Relative x Displacement (mm)", fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(35, 90)
        plt.ylim(-0.1, 1.3)
        # plt.title("Displacement by genotype (score 0)")
        plt.legend(fontsize=16, loc="lower right")

    def plot_stress_by_geno(self):
        """Plot relative effective stress change over time for each genotype."""
        plt.figure()
        for geno in self.genotypes:
            df = self.data.get((geno, 0))
            if df is not None:
                time_years = (df["x"] - 1.1) * 90
                idx = (df["x"] - 1.1).abs().argmin()
                homeostasis_val = df["eff_stress"].iloc[idx]
                stress = (df["eff_stress"] - homeostasis_val) / homeostasis_val
                plt.plot(time_years, stress, label=f"{geno} score 0")
        plt.xlabel("Time (years)", fontsize=20)
        plt.ylabel("Relative Effective Stress Change", fontsize=20)
        plt.xlim(35, 90)
        # plt.ylim(0, 0.5)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        # plt.title("Effective stress by genotype (score 0)")
        plt.legend(fontsize=16, loc="lower right")

    def plot_displacement_by_score(self):
        """Plot relative x displacement over time for each polygenic score."""
        plt.figure()
        for score in self.scores:
            df = self.data.get(("TC", score))
            if df is not None:
                time_years = (df["x"] - 1.1) * 90
                idx = (df["x"] - 1.1).abs().argmin()
                homeostasis_val = df["x_disp"].iloc[idx]
                disp = (df["x_disp"] - homeostasis_val) / homeostasis_val
                plt.plot(time_years, disp, label=f"TC score {score}")
        plt.xlabel("Time (years)", fontsize=20)
        plt.ylabel("Relative x Displacement (mm)", fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(35, 90)
        plt.ylim(-0.1, 1.6)
        # plt.title("Relative displacement by score (TC genotype)")
        plt.legend(fontsize=16, loc="lower right")

    def plot_stress_by_score(self):
        """Plot relative effective stress change over time for each polygenic score."""
        plt.figure()
        for score in self.scores:
            df = self.data.get(("TC", score))
            if df is not None:
                # baseline_stress = df.loc[df["x"] == 1.1, "eff_stress"].values[0]
                idx = (df["x"] - 1.1).abs().argmin()
                homeostasis_val = df["eff_stress"].iloc[idx]
                time_years = (df["x"] - 1.1) * 90
                stress = (df["eff_stress"] - homeostasis_val) / homeostasis_val
                plt.plot(time_years, stress, label=f"TC score {score}")
        plt.xlabel("Time (years)", fontsize=20)
        plt.ylabel("Relative Effective Stress Change", fontsize=20)
        plt.xlim(35, 90)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        # plt.title("TC effective stress by score")
        plt.legend(fontsize=16, loc="lower right")

    def main(self):
        self.plot_displacement_by_geno()
        self.plot_stress_by_geno()
        self.plot_displacement_by_score()
        self.plot_stress_by_score()
        plt.show()


if __name__ == "__main__":
    plotter = FEBioPlotter()
    plotter.main()
