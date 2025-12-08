import numpy as np
from scipy.optimize import fsolve

from aneurysm_sim.model import functions

def simulate_arterial_stress_and_pressure(params):

    n = 235
    sv_stretch_var = 0.55 + 0.01 * np.arange(n) # Stretch variable for simulation

    sv_stress_var_elastin = np.zeros(n)
    sv_stress_var_collagen_me = np.zeros(n)
    sv_stress_var_collagen_ad = np.zeros(n)
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
        sv_stress_var_elastin[i] = functions.sigma_elastin(stretch, params)
        sv_stress_var_collagen[i] = functions.sigma_collagen(stretch, params)
        sv_stress_var_collagen_me[i] = functions.sigma_collagen_me(stretch, params)
        sv_stress_var_collagen_ad[i] = functions.sigma_collagen_ad(stretch, params)
        # sv_stress_var_muscle_a[i] = functions.sigma_muscle_a(stretch, params)
        # sv_stress_var_muscle_p[i] = functions.sigma_muscle_p(stretch, params)
        sv_stress_var_muscle_t[i] = functions.sigma_muscle_total(stretch, params)

        sv_stress_var_total[i] = max((
            sv_stress_var_elastin[i]
            + sv_stress_var_collagen[i]
        ), 0)

        # Pressures
        sv_pressure_var_elastin[i] = max(functions.pressure_elastin(stretch, params), 0)
        sv_pressure_var_collagen[i] = functions.pressure_collagen(stretch, params)
        sv_pressure_var_collagen_me[i] = functions.pressure_collagen_me(stretch, params)
        sv_pressure_var_collagen_ad[i] = functions.pressure_collagen_ad(stretch, params)
        # sv_pressure_var_muscle_p[i] = max(functions.pressure_muscle_p(stretch, params), 0)
        sv_pressure_var_muscle_a[i] = max(functions.pressure_muscle_a(stretch, params), 0)
        sv_pressure_var_muscle[i] = sv_pressure_var_muscle_a[i]

        sv_pressure_var[i] = sv_pressure_var_elastin[i] + sv_pressure_var_collagen_me[i] + sv_pressure_var_collagen_ad[i] + sv_pressure_var_muscle[i]                                            
        
    return {
        # Diameter variable
        'sv_diam_var': sv_diam_var,

        # Stretch variable
        'sv_stretch_var': sv_stretch_var,
        
        # Stress variables
        'sv_stress_var_elastin': sv_stress_var_elastin,
        'sv_stress_var_collagen_me': sv_stress_var_collagen_me,
        'sv_stress_var_collagen_ad': sv_stress_var_collagen_ad,
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

def simulate_aneurysm(params, genotype = None, treatment = False, dt = 0.0069): # dt in years, step independence achieved at 0.0069 from paper
    """
    """
    steps = int(params.t_sim / dt) + 1
    time = np.linspace(0, params.t_sim, steps)

    # Medial collagen, elastin, elastases, collagenases, and immune cells
    elastin_me = np.ones(steps)
    collagen_me = np.ones(steps)
    elastases = np.zeros(steps)
    collagenases = np.zeros(steps)
    immune_cells = np.zeros(steps)

    # Adventitial collagen, fibroblast, procollagen, collagenase, zymogen, TIMP, and TGF-beta                            
    fibroblast = np.zeros(steps)
    muscle_cells = np.zeros(steps)
    collagen_ad = np.zeros(steps)
    procollagen = np.zeros(steps)
    collagenase = np.zeros(steps)
    zymogen = np.zeros(steps)
    timp = np.zeros(steps)
    latent_tgf_beta = np.zeros(steps)
    active_tgf_beta = np.zeros(steps)

    # Mechanical parameters
    diameter = np.zeros(steps)
    lambda_c_max = np.zeros(steps)
    lambda_c_min = np.zeros(steps)
    lambda_c_mode = np.zeros(steps)
    lambda_att_max = np.zeros(steps)
    lambda_att_min = np.zeros(steps)
    lambda_att_mode = np.zeros(steps)
    lambda_rec_max = np.zeros(steps)
    lambda_rec_min = np.zeros(steps)
    lambda_rec_mode = np.zeros(steps)
    lambda_sys_array = np.zeros(steps)

    # Initialize variables
    collagen_me[0] = params.init_collagen_me
    elastin_me[0] = params.init_elastin_me
    collagenases[0] = params.init_collagenases
    elastases[0] = params.init_elastases
    immune_cells[0] = params.i_0 
    fibroblast[0] = params.init_fibroblast
    muscle_cells[0] = params.init_muscle_cells
    collagen_ad[0] = params.init_collagen_ad
    procollagen[0] = params.init_procollagen
    collagenase[0] = params.init_collagenase
    zymogen[0] = params.init_zymogen
    timp[0] = params.init_timp 
    latent_tgf_beta[0] = params.init_latent_tgf_beta 
    active_tgf_beta[0] = params.init_active_tgf_beta 
    tgf_beta_level = functions.get_latent_tgf_beta_level(params, genotype)

    # Initialize collagen attachment and recruitment stretch parameters
    lambda_att_max[0]  = params.c_att_max_ad
    lambda_att_min[0]  = params.c_att_min_ad
    lambda_att_mode[0] = params.c_att_mod_ad
    lambda_rec_max[0]  = params.c_rec_max_ad 
    lambda_rec_min[0]  = params.c_rec_min_ad
    lambda_rec_mode[0] = params.c_rec_mod_ad 
    lambda_c_max[0] = params.c_lambda_sys / lambda_rec_min[0] 
    lambda_c_min[0] = params.c_lambda_sys / lambda_rec_max[0]
    lambda_c_mode[0] = params.c_lambda_sys/  lambda_rec_mode[0]
    diameter[0] = 2*params.c_radius_tzero*params.c_lambda_sys
    lambda_sys_array[0] = params.c_lambda_sys
    lambd_c_max_history = [lambda_c_max[0]] 

    for i in range(1, steps): 
        t = time[i]

        # Calculate the lambda_sys using fsolve to find the root of the force balance equation
        lambda_sys = fsolve(functions.force_balance_equation, [lambda_sys_array[i-1]], args=(elastin_me[i-1], collagen_me[i-1], collagen_ad[i-1], muscle_cells[i-1], params,))[0]
        lambda_sys_array[i] = lambda_sys
        
        lambda_c_max[i] = lambda_sys / lambda_rec_min[i-1]
        lambda_c_min[i] = lambda_sys / lambda_rec_max[i-1]
        lambda_c_mode[i] = lambda_sys / lambda_rec_mode[i-1]
        lambd_c_max_history.append(lambda_c_max[i])

        lambda_att_max[i] = functions.calculate_max_attachment_stretch(lambd_c_max_history, dt, i, params)
        lambda_att_min[i] = functions.calculate_min_attachment_stretch(lambda_att_max[i], params)
        lambda_att_mode[i] = functions.calculate_mode_attachment_stretch(lambda_att_min[i], lambda_att_max[i], params)

        alpha = functions.alpha_rate(fibroblast[i-1], collagen_ad[i-1], collagenase[i-1], params)

        lambda_rec_min[i] = lambda_rec_min[i-1] + dt * functions.d_collagen_min_recruitment_stretch_ad_dt(alpha, lambda_c_max[i], lambda_att_max[i])
        lambda_rec_max[i] = lambda_rec_max[i-1] + dt * functions.d_collagen_max_recruitment_stretch_ad_dt(alpha, lambda_c_min[i], lambda_att_min[i])
        lambda_rec_mode[i] = lambda_rec_mode[i-1] + dt * functions.d_collagen_mode_recruitment_stretch_ad_dt(alpha, lambda_c_mode[i], lambda_att_mode[i])

        diameter[i] = 2 * params.c_radius_tzero * lambda_sys

        if t < params.t_i0:
            collagen_me[i] = collagen_me[0]
            collagen_ad[i] = collagen_ad[0]
            elastin_me[i] = elastin_me[0]
            collagenases[i] = collagenases[0]
            elastases[i] = elastases[0]
            immune_cells[i] = immune_cells[0]
            fibroblast[i] = fibroblast[0]
            muscle_cells[i] = muscle_cells[0]
            procollagen[i] = procollagen[0]
            collagenase[i] = collagenase[0]
            zymogen[i] = zymogen[0]
            timp[i] = timp[0]
            latent_tgf_beta[i] = latent_tgf_beta[0]
            active_tgf_beta[i] = active_tgf_beta[0]
            lambda_sys_array[i] = lambda_sys_array[0]

        else: 
            immune_cells[i] = functions.calculate_immune_cell_level(t, params)
            collagen_me[i] = collagen_me[i-1] + dt * functions.d_medial_collagen_dt(collagenases[i-1], collagen_me[i-1], params)
            collagen_ad[i] = collagen_ad[i-1] + dt * functions.d_collagen_dt(procollagen[i-1], collagenase[i-1], collagen_ad[i-1], params)
            elastin_me[i] = elastin_me[i-1] + dt * functions.d_medial_elastin_dt(elastases[i-1], elastin_me[i-1], params)
            collagenases[i] = collagenases[i-1] + dt * functions.d_collagenases_dt(immune_cells[i-1], collagenases[i-1], params)
            elastases[i] = elastases[i-1] + dt * functions.d_elastases_dt(immune_cells[i-1], elastases[i-1], params)
            muscle_cells[i] = muscle_cells[i-1] + dt * functions.d_muscle_cells_dt(lambda_sys, muscle_cells[i-1], params)
            fibroblast[i] = fibroblast[i-1] + dt * functions.d_fibroblast_dt(active_tgf_beta[i-1], fibroblast[i-1], params)
            procollagen[i] = procollagen[i-1] + dt * functions.d_procollagen_dt(active_tgf_beta[i-1], fibroblast[i-1], procollagen[i-1], params)
            collagenase[i] = collagenase[i-1] + dt * functions.d_collagenase_dt(collagenase[i-1], zymogen[i-1], timp[i-1], params)
            zymogen[i] = zymogen[i-1] + dt * functions.d_zymogen_dt(active_tgf_beta[i-1], fibroblast[i-1], zymogen[i-1], params)
            timp[i] = timp[i-1] + dt * functions.d_timp_dt(active_tgf_beta[i-1], fibroblast[i-1], collagenase[i-1], timp[i-1], params)
            latent_tgf_beta[i] = latent_tgf_beta[i-1] + dt * functions.d_latent_tgf_beta_dt(active_tgf_beta[i-1], latent_tgf_beta[i-1],
                                                                                fibroblast[i-1], collagen_ad[i-1], lambda_c_max[i], lambda_att_max[i-1], tgf_beta_level, params)
            active_tgf_beta[i] = active_tgf_beta[i-1] + dt * functions.d_active_tgf_beta_dt(active_tgf_beta[i-1], latent_tgf_beta[i-1], fibroblast[i-1], lambda_c_max[i], lambda_att_max[i-1], params)

            if treatment and abs(t - params.t_treat) < dt:
                active_tgf_beta[i] += params.tgf_spike_amount

    return {
        'time': time,
        'elastin_me': elastin_me,
        'collagen_me': collagen_me,
        'elastases': elastases,
        'collagenases': collagenases,
        'immune_cells': immune_cells,
        'fibroblast': fibroblast,
        'muscle_cells': muscle_cells,
        'collagen_ad': collagen_ad,
        'procollagen': procollagen,
        'collagenase': collagenase,
        'zymogen': zymogen,
        'timp': timp,
        'latent_tgf_beta': latent_tgf_beta,
        'active_tgf_beta': active_tgf_beta,
        'diameter': diameter,
        'lambda_c_max': lambda_c_max,
        'lambda_c_min': lambda_c_min,
        'lambda_c_mode': lambda_c_mode,
        'lambda_att_max': lambda_att_max,
        'lambda_att_min': lambda_att_min,
        'lambda_att_mode': lambda_att_mode,
        'lambda_rec_max': lambda_rec_max,
        'lambda_rec_min': lambda_rec_min,
        'lambda_rec_mode': lambda_rec_mode,
        'lambd_c_max_history': lambd_c_max_history,
        'lambda_sys': lambda_sys_array,
        'final_lambda_sys': lambda_sys_array[-1],
    }

