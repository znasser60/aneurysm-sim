import numpy as np 

# Vasospasm mechanical equations
def v_lambda_collagen(x, params): 
    return x / params.c_lambda_elastin # TODO: Original code says c_rec_collagen, but constant does not exist. 

def v_lambda_muscle(x, params):
    return x / params.c_rec_muscle # TODO: Should this be c_lambda_muscle?

def v_m(x, params): 
    return x / params.c_rec_muscle

def v_ge_collagen(x, params):
    return (v_lambda_collagen(x, params)**2 - 1) / 2

def v_ge_muscle(x, params):
    return (v_lambda_muscle(x, params)**2 - 1) / 2

def v_ge(x, params):
    return (x**2 - 1) / 2

def v_sigma_elastin(x, params):
    return x**2 * params.c_k_elastin * (1 - (1 / (params.c_lambda_z**2 * x**4)))

def v_sigma_muscle_p(x, params):
    return v_lambda_muscle(x, params)**2 * params.c_k_muscle_p * (1 - 1 / (params.c_lambda_z**2 * v_lambda_muscle(x, params)**4))

def v_sigma_muscle_a(x, params):
    return params.c_vasodil_conc * params.c_k_muscle_a * v_m(x, params) * (1 - ((params.c_musc_mean - v_m(x, params)) / (params.c_musc_mean - params.c_musc_min))**2)

def v_sigma_muscle_t(x, params):
    return v_sigma_muscle_a(x, params) + v_sigma_muscle_p(x, params)

def v_sigma_collagen_me_0(x):
    return 0 * x

def v_sigma_collagen_me_ac(x, params):
    return x * params.v_gamma_me * 2 * ((x + params.v_a_me) * np.log(x / params.v_a_me) + 2 * (params.v_a_me - x))

def v_sigma_collagen_me_cb(x, params):
    term1 = (x + params.v_a_me) * np.log(params.v_c_me / params.v_a_me) + params.v_a_me - params.v_c_me + ((params.v_a_me - params.v_c_me) / params.v_c_me) * x
    term2 = (x + params.v_b_me) * np.log(x / params.v_c_me) + params.v_b_me + params.v_c_me - ((params.v_b_me + params.v_c_me) / params.v_c_me) * x
    return x * params.v_gamma_me * 2 * term1 - x * params.v_delta_me * 2 * term2

def v_sigma_collagen_me_b(x, params):
    term1 = (x + params.v_a_me) * np.log(params.v_c_me / params.v_a_me) + params.v_a_me - params.v_c_me + ((params.v_a_me - params.v_c_me) / params.v_c_me) * x
    term2 = (x + params.v_b_me) * np.log(params.v_b_me / params.v_c_me) - params.v_b_me + params.v_c_me - ((params.v_b_me - params.v_c_me) / params.v_c_me) * x
    return x * params.v_gamma_me * 2 * term1 - x * params.v_delta_me * 2 * term2

def v_sigma_collagen_me(x, params):
    if x < params.v_a_me:
        return v_sigma_collagen_me_0(x)
    elif x < params.v_c_me:
        return v_sigma_collagen_me_ac(x, params)
    elif x <= params.v_b_me:
        return v_sigma_collagen_me_cb(x, params)
    else:
        return v_sigma_collagen_me_b(x, params)

def v_sigma_collagen_ad_0(x):
    return 0 * x

def v_sigma_collagen_ad_ac(x, params):
    return x * params.v_gamma_ad * 2 * ((x + params.v_a_ad) * np.log(x / params.v_a_ad) + 2 * (params.v_a_ad - x))

def v_sigma_collagen_ad_cb(x, params):
    term1 = (x + params.v_a_ad) * np.log(params.v_c_ad / params.v_a_ad) + params.v_a_ad - params.v_c_ad + ((params.v_a_ad - params.v_c_ad) / params.v_c_ad) * x
    term2 = (x + params.v_b_ad) * np.log(x / params.v_c_ad) + params.v_b_ad + params.v_c_ad - ((params.v_b_ad + params.v_c_ad) / params.v_c_ad) * x
    return x * params.v_gamma_ad * 2 * term1 - x * params.v_delta_ad * 2 * term2

def v_sigma_collagen_ad_b(x, params):
    term1 = (x + params.v_a_ad) * np.log(params.v_c_ad / params.v_a_ad) + params.v_a_ad - params.v_c_ad + ((params.v_a_ad - params.v_c_ad) / params.v_c_ad) * x
    term2 = (x + params.v_b_ad) * np.log(params.v_b_ad / params.v_c_ad) - params.v_b_ad + params.v_c_ad - ((params.v_b_ad - params.v_c_ad) / params.v_c_ad) * x
    return x * params.v_gamma_ad * 2 * term1 - x * params.v_delta_ad * 2 * term2

def v_sigma_collagen_ad(x, params): 
    if x < params.v_a_ad:
        return v_sigma_collagen_ad_0(x)
    elif x < params.v_c_ad:
        return v_sigma_collagen_ad_ac(x, params)
    elif x <= params.v_b_ad:
        return v_sigma_collagen_ad_cb(x, params)
    else:
        return v_sigma_collagen_ad_b(x, params)

def v_sigma_collagen(x, params):
    return v_sigma_collagen_me(x, params) + v_sigma_collagen_ad(x, params)

def v_pres_prefactor(x, params): 
    return params.c_thickness_tzero / (params.c_radius_tzero * params.c_lambda_z * x**2)

def v_pressure_ECM(x, params): 
    return v_pres_prefactor(x, params) * (v_sigma_elastin(x, params) + v_sigma_collagen_me(x, params) + v_sigma_collagen_ad(x, params) + v_sigma_muscle_t(x, params))

def v_pressure_EC(x, params):
    return v_pres_prefactor(x, params) * (v_sigma_elastin(x, params) + v_sigma_collagen_me(x, params) + v_sigma_collagen_ad(x, params))

def v_pressure_EM(x, params):
    return v_pres_prefactor(x, params) * (v_sigma_elastin(x, params) + v_sigma_muscle_t(x, params))

def v_pressure_E(x, params):
    return v_pres_prefactor(x, params) * v_sigma_elastin(x, params)

def v_pressure_elastin(x, params):
    return v_pres_prefactor(x, params) * v_sigma_elastin(x, params)

def v_pressure_collagen(x, params):
    return v_pres_prefactor(x, params) * v_sigma_collagen(x, params)

def v_pressure_muscle(x, params):
    return v_pres_prefactor(x, params) * v_sigma_muscle_t(x, params)

def v_pressure_muscle_a(x, params):
    return v_pres_prefactor(x, params) * v_sigma_muscle_a(x, params)

def v_pressure_muscle_p(x, params):
    return v_pres_prefactor(x, params) * v_sigma_muscle_p(x, params)

def v_pressure_collagen_me(x, params):
    return v_pres_prefactor(x, params) * v_sigma_collagen_me(x, params)

def v_pressure_collagen_ad(x, params):
    return v_pres_prefactor(x, params) * v_sigma_collagen_ad(x, params)

# New equations from Aparicio et al. start here: 
# Force balance equation 
def force_balance_equation(lambda_sys_guess, mE_M, mC_M, mC_A, params):             
    """
    Force balance equation for transmural pressure.
    Returns difference between calculated and target pressure.
    """
    lambda_sys = lambda_sys_guess[0]

    # Individual stresses
    stress_elastin = v_sigma_elastin(lambda_sys, params)
    stress_collagen_me = v_sigma_collagen_me(lambda_sys, params)
    stress_collagen_ad = v_sigma_collagen_ad(lambda_sys, params)
    # stress_muscle = v_sigma_muscle_t(lambda_sys, params)
    mass_density_mult = ((mE_M * stress_elastin) + (mC_M * stress_collagen_me) + (mC_A * stress_collagen_ad)) 
    calculated_pressure = (params.c_thickness_tzero / (params.c_radius_tzero * lambda_sys * params.c_lambda_z)) * mass_density_mult

    # Return residual from target pressure
    return calculated_pressure - params.c_pressure_sys 

# Medial degeneration (by immune cells)
def calculate_immune_cell_level(t, params): 
    """
    Immune cell level as a function of time. Equation 7 from Aparicio et al. 2016.
    """
    if t <= params.t_i0: 
        return params.i_0
    else: 
        return params.i_0 + params.i_max * ((t - params.t_i0) / (params.k_i + (t - params.t_i0)))

def d_medial_elastin_dt(elastases, elastin_me, params): 
    """
    Medial elastin ODE: Equation 8 
    """
    return -params.r_e * elastases * elastin_me

def d_medial_collagen_dt(collagenases, collagen_me, params):
    """
    Medial collagen ODE: Equation 8 from Aparicio et al. 2016.
    """
    return -params.r_cm * collagenases * collagen_me

def d_collagenases_dt(immune_cells, collagenases, params): 
    """
    Collagenases ODE: Equation 9 from Aparicio et al. 2016
    Immune cells will come from calculate_immune_cell_level(t, params)
    """
    return params.r_pc1 * immune_cells - params.r_pc2 * collagenases

def d_elastases_dt(immune_cells, elastases, params):
    """
    Elastases ODE: Equation 9 from Aparicio et al. 2016.
    """
    return params.r_pc1 * immune_cells - params.r_pc2 * elastases

# Adventitial collagen growth and remodeling
def f_lambda_fibroblast(lambda_c_max, lambda_att_max):
    """
    Increased stretch of fibroblast cells above homeostatic values leads ot increased production of latent TGF-beta. This models that.
    """
    if lambda_c_max <= lambda_att_max:
        return 0
    return (lambda_c_max - lambda_att_max) / lambda_att_max

def d_fibroblast_dt(tgf_beta, fibroblast, params): 
    """
    Fibroblast ODE: Equation 10 from Apricio et al. 2016.
    Fibroblasts are the cells that produce collagen and elastin in the ECM.
    """
    return (params.r_f1 + params.r_f2 * tgf_beta) * fibroblast - params.r_f3 * fibroblast 

def d_procollagen_dt(tgf_beta, fibroblast, procollagen, params):
    """
    Procollagen ODE: Equation 11 from Apricio et al. 2016.
    Procollagen is the precursor to collagen, which is activated by collagenase.
    """
    return (params.r_p1 + params.r_p2 * tgf_beta) * fibroblast - params.r_p3 * procollagen

def d_collagen_dt(procollagen, collagenase, collagen, params):
    """
    Collagen ODE: Equation 12 from Apricio et al. 2016.
    Collagen is produced from procollagen and is the main structural protein in the ECM.
    """
    return params.r_c1 * procollagen - params.r_c2 * collagenase * collagen

def d_zymogen_dt(tgf_beta, fibroblast, zymogen, params):
    """
    Zymogen ODE: Equation 13 from Apricio et al. 2016.
    Zymogen is the inactive form of collagenase.
    """
    return (params.r_z1 / (1 + params.r_z2 * tgf_beta)) * fibroblast - params.r_z3 * zymogen

def d_collagenase_dt(collagenase, zymogen, timp, params):
    """
    Collagenase ODE: Equation 14 from Apricio et al. 2016.
    Collagenase is the enzyme that activates procollagen to collagen.
    """
    return params.r_ca1 * zymogen - (params.r_ca2 + params.r_ca3 * timp) * collagenase

def d_timp_dt(tgf_beta, fibroblast, collagenase, timp, params):
    """
    TIMP ODE: Equation 15 from Apricio et al. 2016.
    TIMP is the tissue inhibitor of metalloproteinases, which inhibits collagenase.
    """
    return (params.r_i1 + params.r_i2 * tgf_beta) * fibroblast - (params.r_i3 + params.r_i4 * collagenase) * timp

def d_latent_tgf_beta_dt(tgf_beta, latent_tgf_beta, fibroblast, collagen, lambda_c_max, lambda_att_max, tgf_beta_level, params):
    """
    Latent TGF-beta ODE: Equation 16 from Apricio et al. 2016.
    Latent TGF-beta is the inactive form of TGF-beta, which is activated by collagenase.
    """
    term1 = params.r_betal1 * tgf_beta + params.r_betal2 * f_lambda_fibroblast(lambda_c_max, lambda_att_max)
    term2 = 1 + params.r_betal3 * collagen
    term3 = params.r_betal4 + params.r_betal5 * f_lambda_fibroblast(lambda_c_max, lambda_att_max) * fibroblast
    return (term1 / term2) * fibroblast * tgf_beta_level - term3 * latent_tgf_beta       

def d_active_tgf_beta_dt(tgf_beta, latent_tgf_beta, fibroblast, lambda_c_max, lambda_att_max, params):
    """
    Active TGF-beta ODE: Equation 17 from Apricio et al. 2016.
    Active TGF-beta is the active form of TGF-beta, which stimulates fibroblast proliferation and collagen production.
    """
    return (params.r_beta1  + params.r_beta2 * f_lambda_fibroblast(lambda_c_max, lambda_att_max) * fibroblast) * latent_tgf_beta - params.r_beta3 * tgf_beta

# Collagen remodeling I: 
def calculate_max_attachment_stretch(lambda_c_max_history, dt, t_idx, params):
    """
    Numerical implementation of Equation 18 from Aparício et al. (2016)
    - Uses a moving average over T_at years (N steps)
    """
    N = int(params.remodel_time / dt)
    if t_idx < N:
        avg_stretch =  np.mean(lambda_c_max_history[:t_idx+1])
        scale = dt / params.remodel_time
        lambda_att_max = scale * len(lambda_c_max_history[:t_idx+1]) * avg_stretch
    else: 
        window = lambda_c_max_history[t_idx - N + 1 : t_idx + 1]
        scale = dt / params.remodel_time
        lambda_att_max = scale * N * np.mean(window)
    return lambda_att_max

def calculate_min_attachment_stretch(lambda_att_max, params): 
    """
    Eqiattion 19 from Apricio et al. 2016.
    max attachment stretch at time t minus w(t), where w(t) is the width of the attachment stretch distribution at time t.
    For simplicity, we take w(t) = 0.1
    """
    lambda_att_min = lambda_att_max - params.width_att_dist
    return lambda_att_min

def calculate_mode_attachment_stretch(lambda_att_min, lambda_att_max, params): 
    """
    Equation 20 from Apricio et al. 2016.
    min attachment stretch at time t plus s(t)*(min attachment stretch at time t - max attachment stretch at time t), 
    where s(t) is the skew of the distribution at time t. 
    For simplicity, we take s(t) = 0.5
    """
    lambda_att_mode = lambda_att_min + params.skew_att_dist * (lambda_att_min - lambda_att_max)
    return lambda_att_mode

# Collagen remodeling II: 
def alpha_rate(fibroblast, collagen, collagenase, params): 
    """
    Equation 24 from Apricio et al. 2016.
    """
    alpha = params.alpha_init * (fibroblast / collagen) * np.sqrt(collagen * collagenase)
    return alpha

def d_collagen_min_recruitment_stretch_ad_dt(alpha, lambda_c_max, lambda_att_max):
    """
    Equation 21
    """
    return alpha * (lambda_c_max - lambda_att_max) / lambda_att_max

def d_collagen_max_recruitment_stretch_ad_dt(alpha, lambda_c_min, lambda_att_min):
    """
    Equation 22
    """
    return alpha * (lambda_c_min - lambda_att_min) / lambda_att_min

def d_collagen_mode_recruitment_stretch_ad_dt(alpha, lambda_c_mode, lambda_att_mode):
    """
    Equation 23
    """
    return alpha * (lambda_c_mode - lambda_att_mode) / lambda_att_mode

# Add TGF-beta1 protein level function
def get_latent_tgf_beta_level(params, genotype = None):
    levels = params.tgf_beta_levels
    if genotype is None:
        return 1.0
    return levels.get(genotype.upper(), 1.0)

# def get_active_tgf_beta_level(t, params, tgf_beta, treatment = False): 
#     if treatment and t >= params.t_treat: 
#         tgf_beta += params.tgf_spike_amount
#     return tgf_beta
