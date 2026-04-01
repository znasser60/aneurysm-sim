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

def fit_collagen_me_fiber_params_hgo(): 
    """ did not work """
    params = ArterialParameters()
    lam0 = 1.0

    stretches = np.linspace(lam0, 2.0, 100)
    target_stress = np.array(
        [functions.sigma_collagen_me(lam, params) for lam in stretches]
        )
    hgo_stress = lambda lam, k1, k2: (
        4 * k1 * (lam**2) * (lam**2 - 1) * np.exp(k2 * (lam**2 - 1)**2)
    )

    initial_guesses = [1.0,1.0]
    popt, pcov = opt.curve_fit(
        hgo_stress, 
        stretches,
        target_stress, 
        p0=initial_guesses
    )
    
    fitted_k1 = popt[0]
    fitted_k2 = popt[1]
    print(f"Fitted k1: {fitted_k1:.4f}")
    print(f"Fitted k2: {fitted_k2:.4f}")
    fitted_stress = hgo_stress(stretches, fitted_k1, fitted_k2)
    r_squared = 1 - np.sum((target_stress - fitted_stress) ** 2) / np.sum((target_stress - np.mean(target_stress)) ** 2)
    print(f"R-squared: {r_squared:.4f}")
    plt.figure(figsize=(9, 6))
    plt.plot(stretches, target_stress, 'b-', linewidth=3, label='0D Model (Piecewise Logarithmic)')
    plt.plot(stretches, fitted_stress, 'r--', linewidth=2, label=f'HGO Fit (k1={fitted_k1:.2f}, k2={fitted_k2:.2f})')
    plt.axvline(lam0, color='gray', linestyle=':', label=f'Recruitment lam0 ({lam0:.2f})')
    plt.title('Medial Collagen Cauchy Stress Fit (HGO)', fontsize=14, weight='bold')
    plt.xlabel('Circumferential Stretch ($\lambda$)', fontsize=12)
    plt.ylabel('Cauchy Stress (Pa)', fontsize=12)
    plt.legend(fontsize=11)
    plt.show()

def fit_media_stress_hgo_params(): 
    params = ArterialParameters()
    lam0 = 1.0
    stretches = np.linspace(lam0,3, 100)
    hgo_stress = lambda lam, c, k1, k2: (
        c*lam**2 + 4*k1*(lam**2-1)*np.exp(k2*(lam**2-1)**2)
    )
    target_stress = np.array(
        [functions.sigma_collagen_me(lam, params)+functions.sigma_elastin(lam,params)+functions.sigma_muscle_p(lam,params) for lam in stretches]
        )
    
    initial_guesses = [1.0,1.0,1.0]
    lower_limits = [1.0,0.0,0.0]
    upper_limits=[np.inf,np.inf,np.inf]
    popt, pcov = opt.curve_fit(
        hgo_stress, 
        stretches,
        target_stress, 
        p0=initial_guesses, 
        bounds=(lower_limits,upper_limits)
    )
    fitted_c=popt[0]
    fitted_k1=popt[1]
    fitted_k2=popt[2]
    fitted_stress = hgo_stress(stretches, fitted_c, fitted_k1, fitted_k2)
    r_squared = 1 - np.sum((target_stress - fitted_stress) ** 2) / np.sum((target_stress - np.mean(target_stress)) ** 2)
    print(f"Fitted c: {fitted_c:.4f}")
    print(f"Fitted k1: {fitted_k1:.4f}")
    print(f"Fitted k2: {fitted_k2:.4f}")
    print(f"R-squared: {r_squared:.4f}")
    plt.figure(figsize=(9, 6))
    plt.plot(stretches, target_stress, 'b-', linewidth=3, label='0D Model')
    plt.plot(stretches, fitted_stress, 'r--', linewidth=2, label=f'HGO Fit (c={fitted_c:.2f}, k1={fitted_k1:.2f}, k2={fitted_k2:.2f})')
    plt.axvline(lam0, color='gray', linestyle=':', label=f'Recruitment lam0 ({lam0:.2f})')
    plt.title('Cauchy Stress Fit (HGO)', fontsize=14, weight='bold')
    plt.xlabel('Circumferential Stretch ($\lambda$)', fontsize=12)
    plt.ylabel('Cauchy Stress (Pa)', fontsize=12)
    plt.legend(fontsize=11)
    plt.show()

def fit_adventitia_stress_hgo_params():
    params = ArterialParameters()
    lam0 = 1.0
    stretches = np.linspace(lam0,3, 100)
    hgo_stress = lambda lam, c, k1, k2: (
        c*lam**2 + 4*k1*(lam**2-1)*np.exp(k2*(lam**2-1)**2)
    )
    target_stress = np.array(
        [functions.sigma_collagen_ad(lam, params) for lam in stretches]
        )
    
    initial_guesses = [1.0,1.0,1.0]
    lower_limits = [1.0,0.0,0.0]
    upper_limits=[np.inf,np.inf,np.inf]
    popt, pcov = opt.curve_fit(
        hgo_stress, 
        stretches,
        target_stress, 
        p0=initial_guesses, 
        bounds=(lower_limits,upper_limits)
    )
    fitted_c=popt[0]
    fitted_k1=popt[1]
    fitted_k2=popt[2]
    fitted_stress = hgo_stress(stretches, fitted_c, fitted_k1, fitted_k2)
    r_squared = 1 - np.sum((target_stress - fitted_stress) ** 2) / np.sum((target_stress - np.mean(target_stress)) ** 2)
    print(f"Fitted c: {fitted_c:.4f}")
    print(f"Fitted k1: {fitted_k1:.4f}")
    print(f"Fitted k2: {fitted_k2:.4f}")
    print(f"R-squared: {r_squared:.4f}")
    plt.figure(figsize=(9, 6))
    plt.plot(stretches, target_stress, 'b-', linewidth=3, label='0D Model')
    plt.plot(stretches, fitted_stress, 'r--', linewidth=2, label=f'HGO Fit (c={fitted_c:.2f}, k1={fitted_k1:.2f}, k2={fitted_k2:.2f})')
    plt.axvline(lam0, color='gray', linestyle=':', label=f'Recruitment lam0 ({lam0:.2f})')
    plt.title('Ad Collagen Cauchy Stress Fit (HGO)', fontsize=14, weight='bold')
    plt.xlabel('Circumferential Stretch ($\lambda$)', fontsize=12)
    plt.ylabel('Cauchy Stress (Pa)', fontsize=12)
    plt.legend(fontsize=11)
    plt.show()


def fit_collagen_ad_fiber_params():
    params = ArterialParameters()
    lam0 = params.c_rec_mod_ad

    stretches = np.linspace(lam0, 2.0, 100)
    target_stress = np.array([functions.sigma_collagen_ad(lam, params) for lam in stretches])

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
    
    plt.title('Adventitial Collagen Cauchy Stress Fit', fontsize=14, weight='bold')
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

    stretches = np.linspace(1.0, 2.5, 100)
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
    # fit_collagen_me_fiber_parameters()
    # fit_smc_parameters()
    # fit_collagen_me_fiber_params_hgo()
    fit_media_stress_hgo_params()
    fit_adventitia_stress_hgo_params()
    # fit_collagen_ad_fiber_params()