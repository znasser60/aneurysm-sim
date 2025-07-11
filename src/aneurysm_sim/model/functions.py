import numpy as np 

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
    
# def v_sigma_collagen_me(x, params):
#     if x < params.v_a_me:
#         return 0.0
#     elif x < params.v_c_me:
#         return x * params.v_gamma_me * 2 * ((x + params.v_a_me) * np.log(x / params.v_a_me) + 2 * (params.v_a_me - x))
#     elif x <= params.v_b_me:
#         return x * params.v_gamma_me * 2 * ((x + params.v_a_me) * np.log(params.v_c_me / params.v_a_me) + params.v_a_me - params.v_c_me + ((params.v_a_me - params.v_c_me) / params.v_c_me) * x) \
#                - x * params.v_delta_me * 2 * ((x + params.v_b_me) * np.log(x / params.v_c_me) + params.v_b_me + params.v_c_me - ((params.v_b_me + params.v_c_me) / params.v_c_me) * x)
#     else:
#         return x * params.v_gamma_me * 2 * ((x + params.v_a_me) * np.log(params.v_c_me / params.v_a_me) + params.v_a_me - params.v_c_me + ((params.v_a_me -params.v_c_me) / params.v_c_me) * x) \
#                - x * params.v_delta_me * 2 * ((x + params.v_b_me) * np.log(params.v_b_me / params.v_c_me) - params.v_b_me + params.v_c_me - ((params.v_b_me -params.v_c_me) /params.v_c_me) * x)


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

# def v_sigma_collagen_ad(x, params):
#     if x < params.v_a_ad:
#         return 0.0
#     elif x < params.v_c_ad:
#         return x * params.v_gamma_ad * 2 * ((x + params.v_a_ad) * np.log(x / params.v_a_ad) + 2 * (params.v_a_ad - x))
#     elif x <= params.v_b_ad:
#         term1 = (x + params.v_a_ad) * np.log(params.v_c_ad / params.v_a_ad) + params.v_a_ad - params.v_c_ad + ((params.v_a_ad - params.v_c_ad) / params.v_c_ad) * x
#         term2 = (x + params.v_b_ad) * np.log(x / params.v_c_ad) + params.v_b_ad + params.v_c_ad - ((params.v_b_ad + params.v_c_ad) / params.v_c_ad) * x
#         return x * params.v_gamma_ad * 2 * term1 - x * params.v_delta_ad * 2 * term2
#     else:
#         term1 = (x + params.v_a_ad) * np.log(params.v_c_ad / params.v_a_ad) + params.v_a_ad - params.v_c_ad + ((params.v_a_ad - params.v_c_ad) / params.v_c_ad) * x
#         term2 = (x + params.v_b_ad) * np.log(params.v_b_ad / params.v_c_ad) - params.v_b_ad + params.v_c_ad - ((params.v_b_ad - params.v_c_ad) / params.v_c_ad) * x
#         return x * params.v_gamma_ad * 2 * term1 - x * params.v_delta_ad * 2 * term2

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

# New equations from Aparicio et al. start here 
def calculate_immune_cell_level(t, params): 
    """
    """
    

def f_lambda_fibroblast(lambda_c_max, lambda_att_max):
    """
    Increased stretch of fibroblast cells above homeostatic values leads ot increased production of latent TGF-beta. This models that.
    """
    return (lambda_c_max - lambda_att_max) / lambda_att_max

def d_fibroblast_dt(t, tgf_beta, fibroblast, params): 
    """
    Fibroblast ODE: Equation 10 from Apricio et al. 2016.
    Fibroblasts are the cells that produce collagen and elastin in the ECM.
    """
    return (params.r_f1 + params.r_f2 * tgf_beta) * fibroblast - params.r_f3 * fibroblast 

def d_procollagen_dt(t, tgf_beta, fibroblast, procollagen, params):
    """
    Procollagen ODE: Equation 11 from Apricio et al. 2016.
    Procollagen is the precursor to collagen, which is activated by collagenase.
    """
    return (params.r_p1 + params.r_p2 * tgf_beta) * fibroblast - params.r_p3 * procollagen

def d_collagen_dt(t, procollagen, collagenase, collagen, params):
    """
    Collagen ODE: Equation 12 from Apricio et al. 2016.
    Collagen is produced from procollagen and is the main structural protein in the ECM.
    """
    return params.r_c1 * procollagen - params.r_c2 * collagenase * collagen

def d_zymogen_dt(t, tgf_beta, fibroblast, zymogen, params):
    """
    Zymogen ODE: Equation 13 from Apricio et al. 2016.
    Zymogen is the inactive form of collagenase.
    """
    return (params.r_z1 / (1 + params.r_z2 * tgf_beta)) * fibroblast - params.r_z3 * zymogen

def d_collagenase_dt(t, collagenase, zymogen, timp, params):
    """
    Collagenase ODE: Equation 14 from Apricio et al. 2016.
    Collagenase is the enzyme that activates procollagen to collagen.
    """
    return params.r_ca1 * zymogen - (params.r_ca2 + params.r_ca3 * timp) * collagenase

def d_timp_dt(t, tgf_beta, fibroblast, collagenase, timp, params):
    """
    TIMP ODE: Equation 15 from Apricio et al. 2016.
    TIMP is the tissue inhibitor of metalloproteinases, which inhibits collagenase.
    """
    return (params.r_l1 + params.r_l2 * tgf_beta) * fibroblast - (params.r_l3 + params.r_l4 * collagenase) * timp

def d_latent_tgf_beta_dt(t, tgf_beta, latent_tgf_beta, fibroblast, collagenase, lambda_c_max, lambda_att_max, params):
    """
    Latent TGF-beta ODE: Equation 16 from Apricio et al. 2016.
    Latent TGF-beta is the inactive form of TGF-beta, which is activated by collagenase.
    """
    term1 = params.r_betal1 * tgf_beta + params.r_betal2 * f_lambda_fibroblast(lambda_c_max, lambda_att_max)
    term2 = 1 + params.r_betal3 * collagenase
    term3 = params.r_betal4 + params.r_betal5 * f_lambda_fibroblast(lambda_c_max, lambda_att_max) * fibroblast
    return (term1 / term2) * fibroblast - term3 * latent_tgf_beta       

def d_active_tgf_beta_dt(t, tgf_beta, latent_tgf_beta, fibroblast, lambda_c_max, lambda_att_max, params):
    """
    Active TGF-beta ODE: Equation 17 from Apricio et al. 2016.
    Active TGF-beta is the active form of TGF-beta, which stimulates fibroblast proliferation and collagen production.
    """
    return (params.r_beta1  + params.r_beta2 * f_lambda_fibroblast(lambda_c_max, lambda_att_max) * fibroblast) * latent_tgf_beta - params.r_beta3 * tgf_beta

def calculate_max_attachment_stretch(lambda_c_max_history, dt, t_idx, params):
    """
    Numerical implementation of Equation 18 from Aparício et al. (2016)
    - Uses a moving average over T_at years (N steps)
    """
    N = int(params.remodel_time / dt)
    if t_idx < N:
        return np.mean(lambda_c_max_history[:t_idx+1])
    else:
        return np.mean(lambda_c_max_history[t_idx - N + 1:t_idx + 1])

def calculate_min_attachment_stretch(lambda_att_max, params): 
    """
    Eqiattion 19 from Apricio et al. 2016.
    max attachment stretch at time t minus w(t), where w(t) is the width of the attachment stretch distribution at time t.
    For simplicity, we take w(t) = 0.1
    """
    return lambda_att_max - params.width_att_dist

def calculate_mode_attachment_stretch(lambda_att_min, lambda_att_max, params): 
    """
    Equation 20 from Apricio et al. 2016.
    min attachment stretch at time t plus s(t)*(min attachment stretch at time t - max attachment stretch at time t), 
    where s(t) is the skew of the distribution at time t. 
    For simplicity, we take s(t) = 0.5
    """
    return lambda_att_min + params.skew_att_dist * (lambda_att_min - lambda_att_max)

# Collagen remodeling II: 
def alpha_rate(params, fibroblast, collagen, collagenase): 
    """
    Equation 24 from Apricio et al. 2016.
    """
    return params.alpha_init * (fibroblast / collagen) * np.sqrt(collagen * collagenase)

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

def get_tgf_beta1_protein_level(genotype, treatment = None):
    levels = {"TT": 0.713, "TC": 0.916, "CC": 1.119}
    return levels.get(genotype.upper(), 1.0)









