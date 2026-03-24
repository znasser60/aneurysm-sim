import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import functions

def fit_collagen_me_fiber_parameters(): 
    """Fit the collagen fiber parameters (alpha, ksi, beta) to the 1D model data."""
    
    params = ArterialParameters()
    lam0 = params.c_rec_mod_me

    stretches = np.linspace(lam0, 2.0, 100)
    target_stress = np.array([functions.sigma_collagen_me(lam, params) for lam in stretches])

    fiber_exp_pow_stress = lambda lam, alpha, ksi, beta: (
        2 * ksi * (lam**2) * (lam**2 - lam0**2)**(beta - 1) * np.exp(alpha * (lam**2 - lam0**2)**beta)
    )

    initial_guesses = [10.0, params.c_k_collagen, 1.6]
    lower_bounds = [0.0, 0.0, 1.0]
    upper_bounds = [np.inf, np.inf, 2.0]

    popt, pcov = opt.curve_fit(
        fiber_exp_pow_stress, 
        stretches, 
        target_stress, 
        p0=initial_guesses, 
        bounds=(lower_bounds, upper_bounds)
    )

    fitted_alpha = popt[0]
    fitted_ksi = popt[1]
    fitted_beta = popt[2]
    print(f"Fitted alpha: {fitted_alpha:.4f}")
    print(f"Fitted ksi: {fitted_ksi:.4f}")
    print(f"Fitted beta: {fitted_beta:.4f}")

    fitted_stress = fiber_exp_pow_stress(stretches, fitted_alpha, fitted_ksi, fitted_beta)
    r_squared = 1 - np.sum((target_stress - fitted_stress) ** 2) / np.sum((target_stress - np.mean(target_stress)) ** 2)
    print(f"R-squared: {r_squared:.4f}")

    plt.figure(figsize=(9, 6))
    plt.plot(stretches, target_stress, 'b-', linewidth=3, label='0D Model (Piecewise Logarithmic)')
    plt.plot(stretches, fitted_stress, 'r--', linewidth=2, label=f'3D FEBio Fit (alpha={fitted_alpha:.2f}, ksi={fitted_ksi:.2f}, beta={fitted_beta:.2f})')
    
    plt.axvline(lam0, color='gray', linestyle=':', label=f'Recruitment lam0 ({lam0:.2f})')
    
    plt.title('Medial Collagen Cauchy Stress Fit', fontsize=14, weight='bold')
    plt.xlabel('Circumferential Stretch ($\lambda$)', fontsize=12)
    plt.ylabel('Cauchy Stress (Pa)', fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.5)
    plt.show()

    return fitted_alpha, fitted_ksi, fitted_beta

def fit_smc_parameters():
    """Fit the smooth muscle cell parameters to the 1D model data."""
    params = ArterialParameters()
    lam0 = params.c_rec_muscle

    stretches = np.linspace(lam0, 2.0, 100)
    target_stress = np.array([functions.sigma_muscle_t(lam, params) for lam in stretches])

    fiber_exp_pow_stress = lambda lam, alpha, ksi, beta: (
        2 * ksi * (lam**2) * (lam**2 - lam0**2)**(beta - 1) * np.exp(alpha * (lam**2 - lam0**2)**beta)
    )

    initial_guesses = [10.0, params.c_k_muscle_a, 2.0]
    lower_bounds = [-np.inf, 0.0, 1.0]
    upper_bounds = [np.inf, np.inf, 2.0]

    popt, pcov = opt.curve_fit(
        fiber_exp_pow_stress, 
        stretches, 
        target_stress, 
        p0=initial_guesses, 
        bounds=(lower_bounds, upper_bounds)
    )

    fitted_alpha = popt[0]
    fitted_ksi = popt[1]
    fitted_beta = popt[2]
    print(f"Fitted alpha: {fitted_alpha:.4f}")
    print(f"Fitted ksi: {fitted_ksi:.4f}")
    print(f"Fitted beta: {fitted_beta:.4f}")

    fitted_stress = fiber_exp_pow_stress(stretches, fitted_alpha, fitted_ksi, fitted_beta)
    r_squared = 1 - np.sum((target_stress - fitted_stress) ** 2) / np.sum((target_stress - np.mean(target_stress)) ** 2)
    print(f"R-squared: {r_squared:.4f}")

    plt.figure(figsize=(9, 6))
    plt.plot(stretches, target_stress, 'b-', linewidth=3, label='0D Model (Piecewise Logarithmic)')
    plt.plot(stretches, fitted_stress, 'r--', linewidth=2, label=f'3D FEBio Fit (alpha={fitted_alpha:.2f}, ksi={fitted_ksi:.2f}, beta={fitted_beta:.2f})')
    
    plt.axvline(lam0, color='gray', linestyle=':', label=f'Recruitment lam0 ({lam0:.2f})')
    
    plt.title('Smooth Muscle Cell Cauchy Stress Fit', fontsize=14, weight='bold')
    plt.xlabel('Circumferential Stretch ($\lambda$)', fontsize=12)
    plt.ylabel('Cauchy Stress (Pa)', fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.5)
    plt.show()

    return fitted_alpha, fitted_ksi, fitted_beta


if __name__ == "__main__":
    fit_collagen_me_fiber_parameters()
    fit_smc_parameters()