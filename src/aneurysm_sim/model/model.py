import numpy as np

from aneurysm_sim.model.functions import (
    v_sigma_elastin, v_sigma_collagen_me, v_sigma_collagen_ad,
    v_sigma_collagen, v_sigma_muscle_a, v_sigma_muscle_p, v_sigma_muscle_t,
    v_pressure_elastin, v_pressure_collagen, v_pressure_collagen_me,
    v_pressure_collagen_ad, v_pressure_muscle_a, v_pressure_muscle_p,
    v_pressure_muscle,
)

from aneurysm_sim.config.parameters import ArterialParameters

def simulate_arterial_stress_and_pressure():
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

    for i in range(n):
        x = sv_stretch_var[i]

        # Stresses
        sv_stress_var_elastin[i] = v_sigma_elastin(x)
        sv_stress_var_collagen[i] = v_sigma_collagen(x)
        sv_stress_var_muscle_a[i] = v_sigma_muscle_a(x)
        sv_stress_var_muscle_p[i] = v_sigma_muscle_p(x)
        sv_stress_var_muscle_t[i] = v_sigma_muscle_t(x)

        sv_stress_var_total[i] = (
            sv_stress_var_elastin[i]
            + sv_stress_var_collagen[i]
            + sv_stress_var_muscle_t[i]
        )

        # Pressures
        sv_pressure_var_elastin[i] = max(v_pressure_elastin(x), 0)
        sv_pressure_var_collagen[i] = max(v_pressure_collagen(x), 0)
        sv_pressure_var_collagen_me[i] = max(v_pressure_collagen_me(x), 0)
        sv_pressure_var_collagen_ad[i] = max(v_pressure_collagen_ad(x), 0)
        sv_pressure_var_muscle_p[i] = max(v_pressure_muscle_p(x), 0)
        sv_pressure_var_muscle_a[i] = max(v_pressure_muscle_a(x), 0)
        sv_pressure_var_muscle[i] = max(sv_pressure_var_muscle_a[i] + sv_pressure_var_muscle_p[i], 0)

        sv_pressure_var[i] = max(
            sv_pressure_var_elastin[i] + sv_pressure_var_collagen_me[i] + sv_pressure_var_collagen_ad[i] + sv_pressure_var_muscle[i], 
            0
            )
        
    # Calculate diameter from stretch variable
    params = ArterialParameters()
    sv_diam_var = 2 * params.c_radius_tzero * 1e3 * sv_stretch_var
        
    return {
        # Diameter variable
        'sv_diam_var': sv_diam_var,
        
        # Stretch variables
        'sv_stretch_var': sv_stretch_var,
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