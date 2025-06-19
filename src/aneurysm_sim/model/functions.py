import numpy as np 

from config import ArterialParameters

params = ArterialParameters()

def v_lambda_collagen(x): 
    return x / params.c_lambda_elastin # TODO: Original code says c_rec_collagen, but constant does not exist. 

def v_lambda_muscle(x):
    return x / params.c_rec_muscle # TODO: Should this be c_lambda_muscle?

def v_m(x): 
    return x / params.c_rec_muscle

def v_ge_collagen(x):
    return (v_lambda_collagen(x)**2 - 1) / 2

def v_ge_muscle(x):
    return (v_lambda_muscle(x)**2 - 1) / 2

def v_ge(x):
    return (x**2 - 1) / 2

def v_sigma_elastin(x):
    return x**2 * params.c_k_elastin * (1 - (1 / (params.c_lambda_z**2 * x**4)))

def v_sigma_muscle_p(x):
    return v_lambda_muscle(x)**2 * params.c_k_muscle_p * (1 - 1 / (params.c_lambda_z**2 * v_lambda_muscle(x)**4))

def v_sigma_muscle_a(x):
    return params.c_vasodil_conc * params.c_k_muscle_a * v_m(x) * (1 - ((params.c_musc_mean - v_m(x)) / (params.c_musc_mean - params.c_musc_min))**2)

def v_sigma_muscle_t(x):
    return v_sigma_muscle_a(x) + v_sigma_muscle_p(x)

def v_sigma_collagen_me_0(x):
    return 0 * x

def v_sigma_collagen_me_ac(x):
    return x * params.v_gamma_me * 2 * ((x + params.v_a_me) * np.log(x / params.v_a_me) + 2 * (params.v_a_me - x))

def v_sigma_collagen_me_cb(x):
    term1 = (x + params.v_a_me) * np.log(params.v_c_me / params.v_a_me) + params.v_a_me - params.v_c_me + ((params.v_a_me - params.v_c_me) / params.v_c_me) * x
    term2 = (x + params.v_b_me) * np.log(x / params.v_c_me) + params.v_b_me + params.v_c_me - ((params.v_b_me + params.v_c_me) / params.v_c_me) * x
    return x * params.v_gamma_me * 2 * term1 - x * params.v_delta_me * 2 * term2

def v_sigma_collagen_me_b(x):
    term1 = (x + params.v_a_me) * np.log(params.v_c_me / params.v_a_me) + params.v_a_me - params.v_c_me + ((params.v_a_me - params.v_c_me) / params.v_c_me) * x
    term2 = (x + params.v_b_me) * np.log(params.v_b_me / params.v_c_me) + params.v_b_me - params.v_c_me + ((params.v_b_me - params.v_c_me) / params.v_c_me) * x
    return x * params.v_gamma_me * 2 * term1 - x * params.v_delta_me * 2 * term2

def v_sigma_collagen_me(x):
    return v_sigma_collagen_me_0(x) * (x < params.v_a_me) + \
            v_sigma_collagen_me_ac(x) * (params.v_a_me <= x) * (x < params.v_c_me) + \
            v_sigma_collagen_me_cb(x) * (params.v_c_me <= x) * (x < params.v_b_me) + \
            v_sigma_collagen_me_b(x) * (x >= params.v_b_me)

def v_sigma_collagen_ad_0(x):
    return 0 * x

def v_sigma_collagen_ad_ac(x):
    return x * params.v_gamma_ad * 2 * ((x + params.v_a_ad) * np.log(x / params.v_a_ad) + 2 * (params.v_a_ad - x))

def v_sigma_collagen_ad_cb(x):
    term1 = (x + params.v_a_ad) * np.log(params.v_c_ad / params.v_a_ad) + params.v_a_ad - params.v_c_ad + ((params.v_a_ad - params.v_c_ad) / params.v_c_ad) * x
    term2 = (x + params.v_b_ad) * np.log(x / params.v_c_ad) + params.v_b_ad + params.v_c_ad - ((params.v_b_ad + params.v_c_ad) / params.v_c_ad) * x
    return x * params.v_gamma_ad * 2 * term1 - x * params.v_delta_ad * 2 * term2

def v_sigma_collagen_ad_b(x):
    term1 = (x + params.v_a_ad) * np.log(params.v_c_ad / params.v_a_ad) + params.v_a_ad - params.v_c_ad + ((params.v_a_ad - params.v_c_ad) / params.v_c_ad) * x
    term2 = (x + params.v_b_ad) * np.log(params.v_b_ad / params.v_c_ad) + params.v_b_ad - params.v_c_ad + ((params.v_b_ad - params.v_c_ad) / params.v_c_ad) * x
    return x * params.v_gamma_ad * 2 * term1 - x * params.v_delta_ad * 2 * term2

def v_sigma_collagen_ad(x):
    return v_sigma_collagen_ad_0(x) * (x < params.v_a_ad) + \
            v_sigma_collagen_ad_ac(x) * (params.v_a_ad <= x) * (x < params.v_c_ad) + \
            v_sigma_collagen_ad_cb(x) * (params.v_c_ad <= x) * (x < params.v_b_ad) + \
            v_sigma_collagen_ad_b(x) * (x >= params.v_b_ad)

def v_sigma_collagen(x): 
    return v_sigma_collagen_me(x) + v_sigma_collagen_ad(x)

def v_pres_prefactor(x): 
    return params.c_thickness_tzero / (params.c_radius_tzero * params.c_lambda_z * x**2)

def v_pressure_ECM(x): 
    return v_pres_prefactor(x) * (v_sigma_elastin(x) + v_sigma_collagen_me(x) + v_sigma_collagen_ad(x) + v_sigma_muscle_t(x))

def v_pressure_EC(x):
    return v_pres_prefactor(x) * (v_sigma_elastin(x) + v_sigma_collagen_me(x) + v_sigma_collagen_ad(x))

def v_pressure_EM(x):
    return v_pres_prefactor(x) * (v_sigma_elastin(x) + v_sigma_muscle_t(x))

def v_pressure_E(x):
    return v_pres_prefactor(x) * v_sigma_elastin(x)

def v_pressure_elastin(x):
    return v_pres_prefactor(x) * v_sigma_elastin(x)

def v_pressure_collagen(x):
    return v_pres_prefactor(x) * v_sigma_collagen(x)

def v_pressure_muscle(x):
    return v_pres_prefactor(x) * v_sigma_muscle_t(x)

def v_pressure_muscle_a(x):
    return v_pres_prefactor(x) * v_sigma_muscle_a(x)

def v_pressure_muscle_p(x):
    return v_pres_prefactor(x) * v_sigma_muscle_p(x)

def v_pressure_collagen_me(x):
    return v_pres_prefactor(x) * v_sigma_collagen_me(x)

def v_pressure_collagen_ad(x):
    return v_pres_prefactor(x) * v_sigma_collagen_ad(x)