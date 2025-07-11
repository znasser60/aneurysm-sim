import numpy as np

from aneurysm_sim.model.functions import (
    v_sigma_elastin, v_sigma_collagen_me, v_sigma_collagen_ad,
    v_sigma_collagen, v_sigma_muscle_a, v_sigma_muscle_p, v_sigma_muscle_t,
    v_pressure_elastin, v_pressure_collagen, v_pressure_collagen_me,
    v_pressure_collagen_ad, v_pressure_muscle_a, v_pressure_muscle_p,
    v_pressure_muscle, get_tgf_beta1_protein_level, calculate_min_attachment_stretch, 
    calculate_mode_attachment_stretch, d_fibroblast_dt, d_active_tgf_beta_dt, 
    d_collagen_dt, d_collagen_max_recruitment_stretch_ad_dt, d_collagen_min_recruitment_stretch_ad_dt, 
    d_collagen_mode_recruitment_stretch_ad_dt, d_collagenase_dt, d_zymogen_dt, d_latent_tgf_beta_dt, 
    d_timp_dt, d_procollagen_dt
)

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
    TEST FUNCTION to simulate elastin degradation and its effect on arterial stress and pressure.
    """
    results_dict = {}
    for degradation_factor in degredation_factors:
        params.c_lambda_elastin *= degradation_factor
        params.c_k_elastin *= degradation_factor
        results = simulate_arterial_stress_and_pressure(params)
        results_dict[f"{int(degradation_factor * 100)}% Elastin"] = results
    return results_dict

def simulate_aneurysm(params, genotype, treatment = None, dt = 0.0069): # dt in years, step independence achieved
    """
    """
    steps = int(params.t_sim / dt) + 1
    time = np.linspace(0, params.t_sim, steps)

    fibroblast = np.zeros(steps)
    collagen = np.zeros(steps)
    procollagen = np.zeros(steps)
    collagenase = np.zeros(steps)
    zymogen = np.zeros(steps)
    timp = np.zeros(steps)
    latent_tgf_beta = np.zeros(steps)
    active_tgf_beta = np.zeros(steps)
    diameter = np.zeros(steps)
    lambda_c_max = np.zeros(steps)
    lambda_att_max = np.zeros(steps)
    lambda_att_min = np.zeros(steps)
    lambda_att_mode = np.zeros(steps)
    lambda_rec_max = np.zeros(steps)
    lambda_rec_min = np.zeros(steps)
    lambda_rec_mode = np.zeros(steps)

    fibroblast[0] = params.init_fibroblast
    collagen[0] = params.init_collagen
    procollagen[0] = params.init_procollagen
    collagenase[0] = params.init_collagenase
    zymogen[0] = params.init_zymogen
    timp[0] = params.init_timp
    latent_tgf_beta[0] = params.init_latent_tgf_beta
    active_tgf_beta[0] = get_tgf_beta1_protein_level(genotype)
    lambda_att_max[0] = 1.0
    lambda_att_min[0] = calculate_min_attachment_stretch(1.0, params)
    lambda_att_mode[0] = calculate_mode_attachment_stretch(lambda_att_min[0], lambda_att_max[0], params)
    lambda_rec_min[0] = lambda_att_min[0]
    lambda_rec_max[0] = lambda_att_max[0]
    lambda_rec_mode[0] = lambda_att_mode[0]
    lambda_c_max[0] = params.c_lambda_elastin / lambda_att_max[0]
    diameter[0] = lambda_att_max[0]

    lambda_c_max_history = [lambda_c_max[0]]

    for i in range(1, steps): 
        # Update TGF-beta levels
        if treatment is not None:
            active_tgf_beta[i] = get_tgf_beta1_protein_level(genotype, treatment)
        else:
            active_tgf_beta[i] = get_tgf_beta1_protein_level(genotype)

        # Update fibroblast, collagen, procollagen, collagenase, zymogen, and TIMP levels using backward Euler method
        fibroblast[i] = fibroblast[i-1] + dt * d_fibroblast_dt(active_tgf_beta[i-1], fibroblast[i-1], params)
        collagen[i] = collagen[i-1] + dt * d_collagen_dt(procollagen[i-1], collagenase[i-1], collagen[i-1], params)
        procollagen[i] = procollagen[i-1] + dt * d_procollagen_dt(active_tgf_beta[i-1], fibroblast[i-1], procollagen[i-1], params)
        collagenase[i] = collagenase[i-1] + dt * d_collagenase_dt(collagenase[i-1], zymogen[i-1], timp[i-1], params)
        zymogen[i] = zymogen[i-1] + dt * d_zymogen_dt(active_tgf_beta[i-1], fibroblast[i-1], zymogen[i-1], params)
        timp[i] = timp[i-1] + dt * d_timp_dt(active_tgf_beta[i-1], fibroblast[i-1], collagenase[i-1], timp[i-1], params)
        latent_tgf_beta[i] = latent_tgf_beta[i-1] + dt * d_latent_tgf_beta_dt(active_tgf_beta[i-1], latent_tgf_beta[i-1],
                                                                              fibroblast[i-1], collagenase[i-1], lambda_c_max, lambda_att_max, params) 
        active_tgf_beta[i] = active_tgf_beta[i-1] + dt * d_active_tgf_beta_dt(active_tgf_beta[i-1], latent_tgf_beta[i-1], fibroblast[i-1], lambda_c_max, lambda_att_max, params) 




































