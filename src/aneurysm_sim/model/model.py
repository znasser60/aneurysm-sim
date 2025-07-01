import numpy as np

from aneurysm_sim.model.functions import (
    v_sigma_elastin, v_sigma_collagen_me, v_sigma_collagen_ad,
    v_sigma_collagen, v_sigma_muscle_a, v_sigma_muscle_p, v_sigma_muscle_t,
    v_pressure_elastin, v_pressure_collagen, v_pressure_collagen_me,
    v_pressure_collagen_ad, v_pressure_muscle_a, v_pressure_muscle_p,
    v_pressure_muscle,
)

from aneurysm_sim.config.parameters import ArterialParameters

def simulate_arterial_stress_and_pressure(params):

    n = 235
    sv_stretch_var = 0.55 + 0.01 * np.arange(n) # Stretch variable for simulation

    sv_stress_var_elastin = np.zeros(n)
    sv_stress_var_collagen = np.zeros(n)
    sv_stress_var_muscle_a = np.zeros(n)
    sv_stress_var_muscle_p = np.zeros(n)
    sv_stress_var_muscle_t = np.zeros(n)
    sv_stress_var_total = np.zeros(n)
    sv_pressure_var = np.zeros(n)
    sv_pressure_var_elastin = np.zeros(n)
    sv_pressure_var_collagen = np.zeros(n)
    sv_pressure_var_muscle = np.zeros(n)
    sv_pressure_var_muscle_a = np.zeros(n)
    sv_pressure_var_muscle_p = np.zeros(n)
    sv_pressure_var_collagen_me = np.zeros(n)
    sv_pressure_var_collagen_ad = np.zeros(n)

    # Calculate diameter from stretch variable
    sv_diam_var = 2 * params.c_radius_tzero * 1e3 * sv_stretch_var

    for i in range(n):
        stretch = sv_stretch_var[i]

        # Stresses
        sv_stress_var_elastin[i] = v_sigma_elastin(stretch, params)
        sv_stress_var_collagen[i] = v_sigma_collagen(stretch, params)
        sv_stress_var_muscle_a[i] = v_sigma_muscle_a(stretch, params)
        sv_stress_var_muscle_p[i] = v_sigma_muscle_p(stretch, params)
        sv_stress_var_muscle_t[i] = v_sigma_muscle_t(stretch, params)

        sv_stress_var_total[i] = (
            sv_stress_var_elastin[i]
            + sv_stress_var_collagen[i]
            + sv_stress_var_muscle_t[i]
        )

        # Pressures
        sv_pressure_var_elastin[i] = max(v_pressure_elastin(stretch, params), 0)
        sv_pressure_var_collagen[i] = v_pressure_collagen(stretch, params)
        sv_pressure_var_collagen_me[i] = v_pressure_collagen_me(stretch, params)
        sv_pressure_var_collagen_ad[i] = v_pressure_collagen_ad(stretch, params)
        sv_pressure_var_muscle_p[i] = max(v_pressure_muscle_p(stretch, params), 0)
        sv_pressure_var_muscle_a[i] = max(v_pressure_muscle_a(stretch, params), 0)
        sv_pressure_var_muscle[i] = sv_pressure_var_muscle_a[i] + sv_pressure_var_muscle_p[i]

        sv_pressure_var[i] = sv_pressure_var_elastin[i] + sv_pressure_var_collagen_me[i] + sv_pressure_var_collagen_ad[i] + sv_pressure_var_muscle[i]

                                                                   
        
    return {
        # Diameter variable
        'sv_diam_var': sv_diam_var,

        # Stretch variable
        'sv_stretch_var': sv_stretch_var,
        
        # Stress variables
        'sv_stress_var_elastin': sv_stress_var_elastin,
        'sv_stress_var_collagen': sv_stress_var_collagen,
        'sv_stress_var_muscle_a': sv_stress_var_muscle_a, 
        'sv_stress_var_muscle_p': sv_stress_var_muscle_p,
        'sv_stress_var_muscle_t': sv_stress_var_muscle_t,
        'sv_stress_var_total': sv_stress_var_total,

        # Pressure variables
        'sv_pressure_var': sv_pressure_var,
        'sv_pressure_var_elastin': sv_pressure_var_elastin,
        'sv_pressure_var_collagen': sv_pressure_var_collagen,
        'sv_pressure_var_muscle': sv_pressure_var_muscle,
        'sv_pressure_var_muscle_a': sv_pressure_var_muscle_a,
        'sv_pressure_var_muscle_p': sv_pressure_var_muscle_p,
        'sv_pressure_var_collagen_me': sv_pressure_var_collagen_me,
        'sv_pressure_var_collagen_ad': sv_pressure_var_collagen_ad,
    }

def simulate_elastin_degredation(params, degredation_factors):
    """
    Test function to simulate elastin degradation and its effect on arterial stress and pressure.
    """
    results_dict = {}
    for degradation_factor in degredation_factors:
        params.c_lambda_elastin *= degradation_factor
        params.c_k_elastin *= degradation_factor
        results = simulate_arterial_stress_and_pressure(params)
        results_dict[f"{int(degradation_factor * 100)}% Elastin"] = results
    return results_dict

def simulate_smooth_muscle_loss(): 
    """
    Test function to simulate smooth muscle loss and its effect on arterial stress and pressure.
    """
    pass



