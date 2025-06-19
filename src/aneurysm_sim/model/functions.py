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

def v_sigma_collagen_me_0(x, params):
    return 0 * x

def v_sigma_collagen_me_ac(x, params):
    return x * params.v_gamma_me * 2 * ((x + params.v_a_me) * np.log(x / params.v_a_me) + 2 * (params.v_a_me - x))

def v_sigma_collagen_me_cb(x, params):
    term1 = (x + params.v_a_me) * np.log(params.v_c_me / params.v_a_me) + params.v_a_me - params.v_c_me + ((params.v_a_me - params.v_c_me) / params.v_c_me) * x
    term2 = (x + params.v_b_me) * np.log(x / params.v_c_me) + params.v_b_me + params.v_c_me - ((params.v_b_me + params.v_c_me) / params.v_c_me) * x
    return x * params.v_gamma_me * 2 * term1 - x * params.v_delta_me * 2 * term2

def v_sigma_collagen_me_b(x, params):
    term1 = (x + params.v_a_me) * np.log(params.v_c_me / params.v_a_me) + params.v_a_me - params.v_c_me + ((params.v_a_me - params.v_c_me) / params.v_c_me) * x
    term2 = (x + params.v_b_me) * np.log(params.v_b_me / params.v_c_me) + params.v_b_me - params.v_c_me + ((params.v_b_me - params.v_c_me) / params.v_c_me) * x
    return x * params.v_gamma_me * 2 * term1 - x * params.v_delta_me * 2 * term2

def v_sigma_collagen_me(x, params):
    return v_sigma_collagen_me_0(x, params) * (x < params.v_a_me) + \
            v_sigma_collagen_me_ac(x, params) * (params.v_a_me <= x) * (x < params.v_c_me) + \
            v_sigma_collagen_me_cb(x, params) * (params.v_c_me <= x) * (x < params.v_b_me) + \
            v_sigma_collagen_me_b(x, params) * (x >= params.v_b_me)

def v_sigma_collagen_ad_0(x, params):
    return 0 * x

def v_sigma_collagen_ad_ac(x, params):
    return x * params.v_gamma_ad * 2 * ((x + params.v_a_ad) * np.log(x / params.v_a_ad) + 2 * (params.v_a_ad - x))

def v_sigma_collagen_ad_cb(x, params):
    term1 = (x + params.v_a_ad) * np.log(params.v_c_ad / params.v_a_ad) + params.v_a_ad - params.v_c_ad + ((params.v_a_ad - params.v_c_ad) / params.v_c_ad) * x
    term2 = (x + params.v_b_ad) * np.log(x / params.v_c_ad) + params.v_b_ad + params.v_c_ad - ((params.v_b_ad + params.v_c_ad) / params.v_c_ad) * x
    return x * params.v_gamma_ad * 2 * term1 - x * params.v_delta_ad * 2 * term2

def v_sigma_collagen_ad_b(x, params):
    term1 = (x + params.v_a_ad) * np.log(params.v_c_ad / params.v_a_ad) + params.v_a_ad - params.v_c_ad + ((params.v_a_ad - params.v_c_ad) / params.v_c_ad) * x
    term2 = (x + params.v_b_ad) * np.log(params.v_b_ad / params.v_c_ad) + params.v_b_ad - params.v_c_ad + ((params.v_b_ad - params.v_c_ad) / params.v_c_ad) * x
    return x * params.v_gamma_ad * 2 * term1 - x * params.v_delta_ad * 2 * term2

def v_sigma_collagen_ad(x, params):
    return v_sigma_collagen_ad_0(x, params) * (x < params.v_a_ad) + \
            v_sigma_collagen_ad_ac(x, params) * (params.v_a_ad <= x) * (x < params.v_c_ad) + \
            v_sigma_collagen_ad_cb(x, params) * (params.v_c_ad <= x) * (x < params.v_b_ad) + \
            v_sigma_collagen_ad_b(x, params) * (x >= params.v_b_ad)

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