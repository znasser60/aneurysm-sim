import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import functions


def fit_media_stress_hgo_params(
    params: ArterialParameters = None, show_plot: bool = False
):
    """
    Fits HGO parameters to the 0D model media Cauchy stress curve.
    Accepts an existing ArterialParameters object to allow patient-specific fitting.
    """
    if params is None:
        params = ArterialParameters()

    lambda_z = params.lambda_z

    kappa_media = 0.27
    gamma_media = 7
    lam0 = 1.0
    stretches = np.linspace(lam0, 1.6, 1000)

    target_stress_Pa = np.array(
        [
            functions.sigma_collagen_me(lam, params)
            + functions.sigma_elastin(lam, params)
            + functions.sigma_muscle_p(lam, params)
            for lam in stretches
        ]
    )
    target_stress_MPa = target_stress_Pa / 1e6

    c_init = (params.k_elastin + params.k_muscle_p) / 2 / 1e6

    def hgo_stress(lam, k1, k2):
        gamma = np.radians(gamma_media)
        kappa = kappa_media

        lambda_r = 1.0 / (lam * lambda_z)

        I1 = lam**2 + lambda_z**2 + lambda_r**2
        I4 = lam**2 * np.cos(gamma) ** 2 + lambda_z**2 * np.sin(gamma) ** 2

        E_val = kappa * (I1 - 3.0) + (1.0 - 3.0 * kappa) * (I4 - 1.0)
        E = np.maximum(E_val, 0.0)

        sigma_matrix = c_init * (lam**2 - lambda_r**2)

        fiber_coeff = (
            kappa * (lam**2 - lambda_r**2)
            + (1.0 - 3.0 * kappa) * (lam**2) * np.cos(gamma) ** 2
        )
        sigma_fiber = 4.0 * k1 * E * np.exp(k2 * E**2) * fiber_coeff

        return sigma_matrix + sigma_fiber

    popt, _ = opt.curve_fit(
        hgo_stress,
        stretches,
        target_stress_MPa,
        p0=[0.1, 0.01],
        bounds=([0, 0], [np.inf, np.inf]),
    )
    k1, k2 = popt

    if show_plot:
        fitted = hgo_stress(stretches, k1, k2)
        r2 = 1 - np.sum((target_stress_MPa - fitted) ** 2) / np.sum(
            (target_stress_MPa - np.mean(target_stress_MPa)) ** 2
        )
        print(
            f"Media — c: {c_init:.6f} MPa, k1: {k1:.6f} MPa, k2: {k2:.6f}, R²: {r2:.4f}"
        )
        print(f"  kappa: {kappa_media}, gamma: {gamma_media}°")

        plt.figure(figsize=(9, 6))
        plt.plot(stretches, target_stress_MPa, "b-", lw=3, label="0D Model (MPa)")
        plt.plot(
            stretches,
            fitted,
            "r--",
            lw=2,
            label=f"HGO Fit (c={c_init:.4f}, k1={k1:.4f}, k2={k2:.4f})",
        )
        plt.title(
            "Media Cauchy Stress — Effective 1D HGO Model", fontsize=14, weight="bold"
        )
        plt.xlabel(r"Circumferential Stretch ($\lambda_\theta$)")
        plt.ylabel("Effective Stress (MPa)")
        plt.legend()
        plt.tight_layout()
        plt.show()

    return c_init, k1, k2


def fit_adventitia_stress_hgo_params(
    params: ArterialParameters = None, show_plot: bool = False
):
    """
    Fits HGO parameters to the 0D model adventitia Cauchy stress curve.
    Accepts an existing ArterialParameters object to allow patient-specific fitting.
    """
    if params is None:
        params = ArterialParameters()

    lambda_z = params.lambda_z
    kappa_ad = 0.32
    gamma_ad = 77
    lam0 = 1.0
    stretches = np.linspace(lam0, 1.6, 1000)

    target_stress_Pa = np.array(
        [functions.sigma_collagen_ad(lam, params) for lam in stretches]
    )
    target_stress_MPa = target_stress_Pa / 1e6
    c_ad = 0.001

    def hgo_stress(lam, k1, k2):
        gamma = np.radians(gamma_ad)
        kappa = kappa_ad

        lambda_r = 1.0 / (lam * lambda_z)
        I1 = lam**2 + lambda_z**2 + lambda_r**2
        I4 = lam**2 * np.cos(gamma) ** 2 + lambda_z**2 * np.sin(gamma) ** 2
        E_val = kappa * (I1 - 3.0) + (1.0 - 3.0 * kappa) * (I4 - 1.0)
        E = np.maximum(E_val, 0.0)
        sigma_matrix = c_ad * (lam**2 - lambda_r**2)
        fiber_coeff = (
            kappa * (lam**2 - lambda_r**2)
            + (1.0 - 3.0 * kappa) * (lam**2) * np.cos(gamma) ** 2
        )
        sigma_fiber = 4.0 * k1 * E * np.exp(k2 * E**2) * fiber_coeff

        return sigma_matrix + sigma_fiber

    popt, _ = opt.curve_fit(
        hgo_stress,
        stretches,
        target_stress_MPa,
        p0=[0.4, 2.0],
        bounds=([0, 0.02], [np.inf, np.inf]),
    )
    k1, k2 = popt

    if show_plot:
        fitted = hgo_stress(stretches, k1, k2)
        r2 = 1 - np.sum((target_stress_MPa - fitted) ** 2) / np.sum(
            (target_stress_MPa - np.mean(target_stress_MPa)) ** 2
        )
        print(
            f"Adventitia — c: {c_ad:.6f} MPa, k1: {k1:.6f} MPa, k2: {k2:.6f}, R²: {r2:.4f}"
        )
        print(f"  kappa: {kappa_ad}, gamma: {gamma_ad}°")

        plt.figure(figsize=(9, 6))
        plt.plot(stretches, target_stress_MPa, "b-", lw=3, label="0D Model (MPa)")
        plt.plot(
            stretches,
            fitted,
            "r--",
            lw=2,
            label=f"HGO Fit (c={c_ad:.4f}, k1={k1:.4f}, k2={k2:.4f})",
        )
        plt.title(
            "Adventitia Cauchy Stress — Effective 1D HGO Model",
            fontsize=14,
            weight="bold",
        )
        plt.xlabel(r"Circumferential Stretch ($\lambda_\theta$)")
        plt.ylabel("Effective Stress (MPa)")
        plt.legend()
        plt.tight_layout()
        plt.show()

    return c_ad, k1, k2


if __name__ == "__main__":
    fit_media_stress_hgo_params(show_plot=True)
    fit_adventitia_stress_hgo_params(show_plot=True)
