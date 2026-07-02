import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "results")
DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data")

class FEBioPlotter:
    # def __init__(self, results_path=RESULTS_PATH):
    #     self.genotypes = ["TT", "TC", "CC"]
    #     self.scores = [0, 4]
    #     self.data = {}

    #     file_name_pattern = re.compile(r"sim_([A-Z]+)_score(\d+)\.log")
    #     for fname in os.listdir(results_path):
    #         match = file_name_pattern.match(fname)
    #         if not match:
    #             continue
    #         geno  = match.group(1)
    #         score = int(match.group(2))
    #         df = self._parse_log(os.path.join(results_path, fname))
    #         if df is not None:
    #             self.data[(geno, score)] = df

    # def _parse_log(self, log_path):
    #     node_rows = []
    #     elem_rows = []
    #     current_type = None
    #     current_time = None

    #     with open(log_path) as f:
    #         for line in f:
    #             line = line.strip()
    #             if not line:
    #                 continue
    #             if line.startswith("Data = ux"):
    #                 current_type = "node"
    #                 continue
    #             if line.startswith("Data = sx"):
    #                 current_type = "elem"
    #                 continue
    #             if line.startswith("Time ="):
    #                 current_time = float(line.split("=")[1].strip())
    #                 continue
    #             if line.startswith("Step") or line.startswith("=") or line.startswith("Data Record"):
    #                 continue
    #             parts = line.split(",")
    #             try:
    #                 vals = [float(x) for x in parts]
    #             except ValueError:
    #                 continue
    #             if current_type == "node" and current_time is not None:
    #                 node_rows.append([current_time] + vals)
    #             elif current_type == "elem" and current_time is not None:
    #                 elem_rows.append([current_time] + vals)

    #     if not node_rows:
    #         return None

    #     node_df = pd.DataFrame(node_rows, columns=["x", "node_id", "ux", "uy", "uz"])
    #     elem_df = pd.DataFrame(elem_rows, columns=["x", "elem_id", "sx", "sy", "sz", "sxy", "syz", "sxz", "s1", "s2", "s3"])

    #     elem_df["eff_stress"] = np.sqrt(0.5 * (
    #         (elem_df["s1"] - elem_df["s2"])**2 +
    #         (elem_df["s2"] - elem_df["s3"])**2 +
    #         (elem_df["s3"] - elem_df["s1"])**2
    #     ))

    #     # drop groupby since only one node/element
    #     node_df = node_df.drop(columns="node_id").rename(columns={"ux": "x_disp"})
    #     elem_df = elem_df[["x", "eff_stress"]]

    #     return node_df.merge(elem_df, on="x")
    
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
                    df = pd.read_csv(file_path, sep=r'\s+').rename(columns={"E1095": meas})
                    dfs.append(df)
                if dfs:
                    merged = dfs[0]
                    for df in dfs[1:]:
                        merged = merged.merge(df, on="x", how="outer")
                    self.data[(geno, score)] = merged.sort_values("x").reset_index(drop=True)

    def plot_displacement_by_geno(self):
        plt.figure()
        for geno in self.genotypes:
            df = self.data.get((geno, 0))
            if df is not None:
                time_years = (df["x"] - 1.1) * 90
                # baseline_disp = df.loc[df["x"] == 1.1, "x_disp"].values[0]
                idx = (df["x"] - 1.1).abs().argmin()
                homeostasis_val = df["x_disp"].iloc[idx]
                disp = (df["x_disp"]-homeostasis_val) / homeostasis_val
                plt.plot(time_years, disp, label=f"{geno} score 0")
        plt.xlabel("Time (years)", fontsize=20)
        plt.ylabel("Relative x Displacement (mm)", fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(35, 90)
        plt.ylim(-0.1, 1.3)
        # plt.title("Displacement by genotype (score 0)")
        plt.legend(fontsize=16, loc='lower right')

    def plot_stress_by_geno(self):
        plt.figure()
        for geno in self.genotypes:
            df = self.data.get((geno, 0))
            if df is not None:
                time_years = (df["x"] - 1.1) * 90
                # baseline_stress = df.loc[df["x"] == 1.1, "eff_stress"].values[0]
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
        plt.legend(fontsize=16, loc='lower right')

    def plot_displacement_by_score(self):
        plt.figure()
        for score in self.scores:
            df = self.data.get(("TC", score))
            if df is not None:
                time_years = (df["x"] - 1.1) * 90
                idx = (df["x"] - 1.1).abs().argmin()
                homeostasis_val = df["x_disp"].iloc[idx]
                disp = (df["x_disp"]-homeostasis_val) / homeostasis_val
                plt.plot(time_years, disp, label=f"TC score {score}")
        plt.xlabel("Time (years)", fontsize=20)
        plt.ylabel("Relative x Displacement (mm)", fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(35, 90)
        plt.ylim(-0.1, 1.6)
        # plt.title("Relative displacement by score (TC genotype)")
        plt.legend(fontsize=16, loc='lower right')

    def plot_stress_by_score(self):
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
        # plt.ylim(0, 0.5)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        # plt.title("TC effective stress by score")
        plt.legend(fontsize=16, loc='lower right')

    def main(self):
        self.plot_displacement_by_geno()
        # self.plot_stress_by_geno()
        self.plot_displacement_by_score()
        self.plot_stress_by_score()
        plt.show()

if __name__ == "__main__":
    plotter = FEBioPlotter()
    plotter.main()