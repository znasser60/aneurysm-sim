import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.optimize import least_squares
from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import functions, model

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
    stretches = np.linspace(lam0,3,1000)
    # hgo_stress = lambda lam, c, k1, k2: (
    #     c*lam**2 + 4*k1*(lam**2-1)*np.exp(k2*(lam**2-1)**2)
    # )
    # hgo_stress = lambda lam, c, k1, k2: (
    #     c*lam**2 + 4*k1*(lam**2-1)*np.exp(k2*(lam**2-1)**2)
    # )
    c = (params.c_k_elastin + params.c_k_muscle_p)/2
    hgo_stress = lambda lam, k1, k2: (
        c * (lam**2) + 4*k1*(lam**2-1)*np.exp(k2*(lam**2-1)**2)
    )
    target_stress = np.array(
        [functions.sigma_collagen_me(lam, params)+functions.sigma_elastin(lam,params)+functions.sigma_muscle_p(lam,params) for lam in stretches]
        )
    
    initial_guesses = [1.0,1.0]
    lower_limits = [0.0, 1e-4]
    upper_limits = [np.inf, np.inf]
    popt, pcov = opt.curve_fit(
        hgo_stress, 
        stretches,
        target_stress, 
        p0=initial_guesses, 
        bounds=(lower_limits,upper_limits)
    )
    fitted_c=c
    fitted_k1=popt[0]
    fitted_k2=popt[1]
    fitted_stress = hgo_stress(stretches, fitted_k1, fitted_k2)
    r_squared = 1 - np.sum((target_stress - fitted_stress) ** 2) / np.sum((target_stress - np.mean(target_stress)) ** 2)
    print(f"Fitted c: {fitted_c}")
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
    stretches = np.linspace(lam0,3, 1000)
    c = 0.001
    hgo_stress = lambda lam, k1, k2: (
        c * (lam**2) + 4*k1*(lam**2-1)*np.exp(k2*(lam**2-1)**2)
    )
    target_stress = np.array(
        [functions.sigma_collagen_ad(lam, params) for lam in stretches]
        )
    
    initial_guesses = [1.0,2.0]
    lower_limits = [0.0, 1e-4]
    upper_limits = [np.inf, np.inf]
    popt, pcov = opt.curve_fit(
        hgo_stress, 
        stretches,
        target_stress, 
        p0=initial_guesses, 
        bounds=(lower_limits,upper_limits)
    )
    fitted_c=c
    fitted_k1=popt[0]
    fitted_k2=popt[1]
    fitted_stress = hgo_stress(stretches, fitted_k1, fitted_k2)
    r_squared = 1 - np.sum((target_stress - fitted_stress) ** 2) / np.sum((target_stress - np.mean(target_stress)) ** 2)
    print(f"Fitted c: {fitted_c}")
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

    stretches = np.linspace(lam0, 2.0, 1000)
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

def calculate_stable_hgo_parameters(c_1D, gamma_deg, lambda_z_max, layer_name):
    """
    Calculates stable FEBio parameters for an HGO material to prevent negative Jacobians.
    
    Parameters:
    c_1D (float): The ground matrix parameter derived from your 1D model (in MPa).
    gamma_deg (float): Fiber angle relative to circumferential direction (degrees).
    lambda_z_max (float): The maximum axial stretch applied before pressure inflation.
    layer_name (str): Name for logging purposes.
    """
    gamma_rad = np.radians(gamma_deg)
    
    # 1. Check Fiber Activation (I4) under pure axial stretch
    # I4 = cos^2(gamma)/lambda_z + lambda_z^2 * sin^2(gamma)
    I4 = (np.cos(gamma_rad)**2 / lambda_z_max) + (lambda_z_max**2 * np.sin(gamma_rad)**2)
    
    fibers_active = I4 > 1.0
    
    print(f"--- Stability Analysis: {layer_name} ---")
    print(f"Initial 1D 'c': {c_1D:.5f} MPa")
    print(f"Fiber I4 at stretch {lambda_z_max}: {I4:.4f}")
    
    # 2. Adjust 'c' for Numerical Stability if Fibers are Slack
    # If fibers go slack during stretch, 'c' alone must support the matrix.
    # An empirical FEA stability threshold for soft tissues is c >= 0.05 MPa.
    c_stable = c_1D
    if not fibers_active:
        print("WARNING: Fibers go SLACK during axial stretch! Ground matrix must support all load.")
        if c_stable < 0.05:
            c_stable = 0.05
            print(f"-> ADJUSTMENT: Boosted 'c' to {c_stable:.5f} MPa to prevent element collapse.")
    else:
        print("Fibers remain engaged. Ground matrix 'c' is safe.")
        
    # 3. Calculate Perfect Bulk Modulus (k)
    # Target a Poisson's ratio of 0.49 to prevent volumetric locking
    nu_target = 0.49
    # k = 4*c * (1+nu) / (3 * (1-2nu))
    k_stable = (4 * c_stable * (1 + nu_target)) / (3 * (1 - 2 * nu_target))
    
    print(f"-> REQUIRED BULK MODULUS <k>: {k_stable:.4f} MPa\n")
    
    return {
        "c_stable": c_stable,
        "k_stable": k_stable,
        "fibers_active": fibers_active
    }
def fit_media_stress_hgo_params2(): 
    params = ArterialParameters()
    
    # 1. Start fitting at 1.0 (no compression!)
    lam0 = 1.0
    stretches = np.linspace(lam0, 3.0, 1000)
    
    # We need the axial stretch to calculate the 3D radial thinning
    lam_z = params.c_lambda_z 
    gamma_deg = 29.0 # Media fiber angle
    gamma_rad = np.radians(gamma_deg)
    
    # Calculate analytical 'c' in Pascals, then convert to MPa for FEBio
    c_Pa = (params.c_k_elastin + params.c_k_muscle_p) / 2 
    c_MPa = c_Pa / 1e6
    
    # 2. The True 3D Continuum Mechanics Equation
    def hgo_3d_stress(lam_theta, k1, k2):
        # Radial thinning (incompressibility)
        lam_r = 1.0 / (lam_theta * lam_z)
        
        # 4th Invariant (Actual stretch on the angled fibers)
        I4 = (lam_theta**2 * np.cos(gamma_rad)**2) + (lam_z**2 * np.sin(gamma_rad)**2)
        
        # Ground matrix stress
        sigma_iso = 2 * c_MPa * (lam_theta**2 - lam_r**2)
        
        # Fiber stress (only if I4 > 1, meaning fibers are in tension)
        # Using np.where to handle arrays seamlessly
        sigma_fib = np.where(
            I4 > 1.0, 
            2 * k1 * (lam_theta**2) * (np.cos(gamma_rad)**2) * (I4 - 1) * np.exp(k2 * (I4 - 1)**2), 
            0.0
        )
        return sigma_iso + sigma_fib

    # Generate target stress in Pa, then convert to MPa
    target_stress_Pa = np.array([
        functions.sigma_collagen_me(lam, params) + 
        functions.sigma_elastin(lam, params) + 
        functions.sigma_muscle_p(lam, params) 
        for lam in stretches
    ])
    target_stress_MPa = target_stress_Pa / 1e6
    
    # Curve Fit
    initial_guesses = [1.0, 1.0]
    lower_limits = [0.0, 1e-4]
    upper_limits = [np.inf, np.inf]
    popt, pcov = opt.curve_fit(
        hgo_3d_stress, 
        stretches,
        target_stress_MPa, 
        p0=initial_guesses, 
        bounds=(lower_limits, upper_limits)
    )
    
    fitted_k1 = popt[0]
    fitted_k2 = popt[1]
    fitted_stress_MPa = hgo_3d_stress(stretches, fitted_k1, fitted_k2)
    
    r_squared = 1 - np.sum((target_stress_MPa - fitted_stress_MPa) ** 2) / np.sum((target_stress_MPa - np.mean(target_stress_MPa)) ** 2)
    
    print("--- MEDIA PARAMETERS FOR FEBIO ---")
    print(f"c  = {c_MPa:.6f} MPa")
    print(f"k1 = {fitted_k1:.6f} MPa")
    print(f"k2 = {fitted_k2:.6f}")
    print(f"R-squared: {r_squared:.4f}")
    
    # Plotting
    plt.figure(figsize=(9, 6))
    plt.plot(stretches, target_stress_MPa, 'b-', linewidth=3, label='0D Model (Target)')
    plt.plot(stretches, fitted_stress_MPa, 'r--', linewidth=2, label=f'HGO 3D Fit (k1={fitted_k1:.3f}, k2={fitted_k2:.3f})')
    plt.title('Media Cauchy Stress Fit (HGO 3D)', fontsize=14, weight='bold')
    plt.xlabel('Circumferential Stretch ($\lambda_\theta$)', fontsize=12)
    plt.ylabel('Cauchy Stress (MPa)', fontsize=12)
    plt.legend(fontsize=11)
    plt.show()

def fit_adventitia_stress_hgo_params2(): 
    params = ArterialParameters()
    
    # 1. Start fitting at 1.0 (no compression!)
    lam0 = 1.0
    stretches = np.linspace(lam0, 3.0, 1000)
    
    # We need the axial stretch to calculate the 3D radial thinning
    lam_z = params.c_lambda_z 
    gamma_deg = 62.0
    gamma_rad = np.radians(gamma_deg)
    
    # Calculate analytical 'c' in Pascals, then convert to MPa for FEBio
    c_Pa = (params.c_k_elastin + params.c_k_muscle_p) / 2 
    c_MPa = 0.001
    
    # 2. The True 3D Continuum Mechanics Equation
    def hgo_3d_stress(lam_theta, k1, k2):
        # Radial thinning (incompressibility)
        lam_r = 1.0 / (lam_theta * lam_z)
        
        # 4th Invariant (Actual stretch on the angled fibers)
        I4 = (lam_theta**2 * np.cos(gamma_rad)**2) + (lam_z**2 * np.sin(gamma_rad)**2)
        
        # Ground matrix stress
        sigma_iso = 2 * c_MPa * (lam_theta**2 - lam_r**2)
        
        # Fiber stress (only if I4 > 1, meaning fibers are in tension)
        # Using np.where to handle arrays seamlessly
        sigma_fib = np.where(
            I4 > 1.0, 
            2 * k1 * (lam_theta**2) * (np.cos(gamma_rad)**2) * (I4 - 1) * np.exp(k2 * (I4 - 1)**2), 
            0.0
        )
        return sigma_iso + sigma_fib

    # Generate target stress in Pa, then convert to MPa
    target_stress_Pa = np.array([
        functions.sigma_collagen_ad(lam, params)
        for lam in stretches
    ])
    target_stress_MPa = target_stress_Pa / 1e6
    
    # Curve Fit
    initial_guesses = [1.0, 1.0]
    lower_limits = [0.0, 1e-4]
    upper_limits = [np.inf, np.inf]
    popt, pcov = opt.curve_fit(
        hgo_3d_stress, 
        stretches,
        target_stress_MPa, 
        p0=initial_guesses, 
        bounds=(lower_limits, upper_limits)
    )
    
    fitted_k1 = popt[0]
    fitted_k2 = popt[1]
    fitted_stress_MPa = hgo_3d_stress(stretches, fitted_k1, fitted_k2)
    
    r_squared = 1 - np.sum((target_stress_MPa - fitted_stress_MPa) ** 2) / np.sum((target_stress_MPa - np.mean(target_stress_MPa)) ** 2)
    
    print("--- MEDIA PARAMETERS FOR FEBIO ---")
    print(f"c  = {c_MPa:.6f} MPa")
    print(f"k1 = {fitted_k1:.6f} MPa")
    print(f"k2 = {fitted_k2:.6f}")
    print(f"R-squared: {r_squared:.4f}")
    
    # Plotting
    plt.figure(figsize=(9, 6))
    plt.plot(stretches, target_stress_MPa, 'b-', linewidth=3, label='0D Model (Target)')
    plt.plot(stretches, fitted_stress_MPa, 'r--', linewidth=2, label=f'HGO 3D Fit (k1={fitted_k1:.3f}, k2={fitted_k2:.3f})')
    plt.title('Media Cauchy Stress Fit (HGO 3D)', fontsize=14, weight='bold')
    plt.xlabel('Circumferential Stretch ($\lambda_\theta$)', fontsize=12)
    plt.ylabel('Cauchy Stress (MPa)', fontsize=12)
    plt.legend(fontsize=11)
    plt.show()

if __name__ == "__main__":
    # fit_collagen_me_fiber_parameters()
    # fit_smc_parameters()
    # fit_collagen_me_fiber_params_hgo()
    fit_media_stress_hgo_params()
    fit_adventitia_stress_hgo_params()
    # fit_media_stress_hgo_params2()
    # fit_adventitia_stress_hgo_params2()
    # fit_collagen_ad_fiber_params()
    # axial_stretch = 1.1

    # Media: c = 0.106 MPa, gamma = 29 deg
    # media_params = calculate_stable_hgo_parameters(
    #     c_1D=0.106, 
    #     gamma_deg=29, 
    #     lambda_z_max=axial_stretch, 
    #     layer_name="Media"
    # )

    # # Adventitia: c = 0.01 MPa, gamma = 62 deg
    # adv_params = calculate_stable_hgo_parameters(
    #     c_1D=0.01, 
    #     gamma_deg=62, 
    #     lambda_z_max=axial_stretch, 
    #     layer_name="Adventitia"
    # )