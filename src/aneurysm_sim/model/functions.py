import numpy as np


# Constituent-specific stretch and stress functions
def lambda_muscle(x, params):
    """Smooth muscle stretch relative to the homeostatic muscle stretch."""
    return x / params.rec_muscle


def lambda_elastin(x, params):
    """Elastin stretch (identity: elastin reference is the wall reference)."""
    return x


def sigma_elastin(x, params):
    """Cauchy stress of elastin (incompressible Neo-hookean)"""
    return (
        lambda_elastin(x, params) ** 2
        * params.k_elastin
        * (1 - (1 / (params.lambda_z**2 * lambda_elastin(x, params) ** 4)))
    )


def sigma_muscle_p(x, params):
    """Cauchy stress of passive smooth muscle (incompressible Neo-hookean)"""
    return (
        lambda_muscle(x, params) ** 2
        * params.k_muscle_p
        * (1 - 1 / (params.lambda_z**2 * lambda_muscle(x, params) ** 4))
    )


def sigma_muscle_a(x, params):
    """Cauchy stress of active smooth muscle response

    Uses the active length-tension relationship where contractile stress
    peaks at musc_mean."""
    return (
        params.vasodil_conc
        * params.k_muscle_a
        * lambda_muscle(x, params)
        * (
            1
            - (
                (params.musc_mean - lambda_muscle(x, params))
                / (params.musc_mean - params.musc_min)
            )
            ** 2
        )
    )


def sigma_muscle_t(x, params):
    """Total Cauchy stress of smooth muscle (passive + active)"""
    return sigma_muscle_a(x, params) + sigma_muscle_p(x, params)


def sigma_collagen_me_0(x):
    """Cauchy stress of medial collagen at zero stretch (no stress)"""
    return 0 * x


def sigma_collagen_me_ac(x, params):
    """Cauchy stress of medial collagen in the ascending section of
    the triangular distribution, between the minimum and mode attachment stretches"""
    return (
        x
        * params.gamma_me
        * 2
        * ((x + params.a_me) * np.log(x / params.a_me) + 2 * (params.a_me - x))
    )


def sigma_collagen_me_cb(x, params):
    """Cauchy stress of medial collagen in the descending section of
    the triangular distribution, between the mode and maximum attachment stretches"""
    term1 = (
        (x + params.a_me) * np.log(params.c_me / params.a_me)
        + params.a_me
        - params.c_me
        + ((params.a_me - params.c_me) / params.c_me) * x
    )
    term2 = (
        (x + params.b_me) * np.log(x / params.c_me)
        + params.b_me
        + params.c_me
        - ((params.b_me + params.c_me) / params.c_me) * x
    )
    return x * params.gamma_me * 2 * term1 - x * params.delta_me * 2 * term2


def sigma_collagen_me_b(x, params):
    """
    Cauchy stress of medial collagen at maximum stretch
    beyond the maximum attachment stretch (all fibers recruited).
    """
    term1 = (
        (x + params.a_me) * np.log(params.c_me / params.a_me)
        + params.a_me
        - params.c_me
        + ((params.a_me - params.c_me) / params.c_me) * x
    )
    term2 = (
        (x + params.b_me) * np.log(params.b_me / params.c_me)
        - params.b_me
        + params.c_me
        - ((params.b_me - params.c_me) / params.c_me) * x
    )
    return x * params.gamma_me * 2 * term1 - x * params.delta_me * 2 * term2


def sigma_collagen_me(x, params):
    """Full piecewise Cauchy stress function of medial collagen
    based on recruitment stages."""
    if x < params.a_me:
        return sigma_collagen_me_0(x)
    elif x < params.c_me:
        return sigma_collagen_me_ac(x, params)
    elif x <= params.b_me:
        return sigma_collagen_me_cb(x, params)
    else:
        return sigma_collagen_me_b(x, params)


def sigma_collagen_ad_0(x):
    """Cauchy stress of adventitial collagen at zero stretch (no stress)"""
    return 0 * x


def sigma_collagen_ad_ac(x, params):
    """Cauchy stress of adventitial collagen in the ascending section of
    the triangular distribution, between the minimum and mode attachment stretches"""
    return (
        x
        * params.gamma_ad
        * 2
        * ((x + params.a_ad) * np.log(x / params.a_ad) + 2 * (params.a_ad - x))
    )


def sigma_collagen_ad_cb(x, params):
    """Cauchy stress of adventitial collagen in the descending section of
    the triangular distribution, between the mode and maximum attachment stretches"""
    term1 = (
        (x + params.a_ad) * np.log(params.c_ad / params.a_ad)
        + params.a_ad
        - params.c_ad
        + ((params.a_ad - params.c_ad) / params.c_ad) * x
    )
    term2 = (
        (x + params.b_ad) * np.log(x / params.c_ad)
        + params.b_ad
        + params.c_ad
        - ((params.b_ad + params.c_ad) / params.c_ad) * x
    )
    return x * params.gamma_ad * 2 * term1 - x * params.delta_ad * 2 * term2


def sigma_collagen_ad_b(x, params):
    """Cauchy stress of adventitial collagen at maximum stretch
    beyond the maximum attachment stretch (all fibers recruited)."""
    term1 = (
        (x + params.a_ad) * np.log(params.c_ad / params.a_ad)
        + params.a_ad
        - params.c_ad
        + ((params.a_ad - params.c_ad) / params.c_ad) * x
    )
    term2 = (
        (x + params.b_ad) * np.log(params.b_ad / params.c_ad)
        - params.b_ad
        + params.c_ad
        - ((params.b_ad - params.c_ad) / params.c_ad) * x
    )
    return x * params.gamma_ad * 2 * term1 - x * params.delta_ad * 2 * term2


def sigma_collagen_ad(x, params):
    """Full piecewise Cauchy stress function of adventitial collagen
    based on recruitment stages."""
    if x < params.a_ad:
        return sigma_collagen_ad_0(x)
    elif x < params.c_ad:
        return sigma_collagen_ad_ac(x, params)
    elif x <= params.b_ad:
        return sigma_collagen_ad_cb(x, params)
    else:
        return sigma_collagen_ad_b(x, params)


def sigma_collagen(x, params):
    """Total collagen Cauchy stress (medial + adventitial)"""
    return sigma_collagen_me(x, params) + sigma_collagen_ad(x, params)


def pres_prefactor(x, params):
    """Laplace prefactor H / (R0 * lambda_z * x^2) converting stress to pressure"""
    return params.thickness_tzero / (params.radius_tzero * params.lambda_z * x**2)


def pressure_ECM(x, params):
    """Total pressure from ECM constituents (elastin, collagen, smooth muscle)"""
    return pres_prefactor(x, params) * (
        sigma_elastin(x, params)
        + sigma_collagen_me(x, params)
        + sigma_collagen_ad(x, params)
        + sigma_muscle_t(x, params)
    )


def pressure_EC(x, params):
    """Total pressure from ECM constituents (elastin, collagen)"""
    return pres_prefactor(x, params) * (
        sigma_elastin(x, params)
        + sigma_collagen_me(x, params)
        + sigma_collagen_ad(x, params)
    )


def pressure_EM(x, params):
    """Total pressure from ECM constituents (elastin, smooth muscle)"""
    return pres_prefactor(x, params) * (
        sigma_elastin(x, params) + sigma_muscle_t(x, params)
    )


def pressure_elastin(x, params):
    """Pressure from elastin constituent"""
    return pres_prefactor(x, params) * sigma_elastin(x, params)


def pressure_collagen(x, params):
    """Pressure from collagen constituent"""
    return pres_prefactor(x, params) * sigma_collagen(x, params)


def pressure_muscle(x, params):
    """Pressure from smooth muscle constituent"""
    return pres_prefactor(x, params) * sigma_muscle_t(x, params)


def pressure_muscle_a(x, params):
    """Pressure from active smooth muscle constituent"""
    return pres_prefactor(x, params) * sigma_muscle_a(x, params)


def pressure_muscle_p(x, params):
    """Pressure from passive smooth muscle constituent"""
    return pres_prefactor(x, params) * sigma_muscle_p(x, params)


def pressure_collagen_me(x, params):
    """Pressure from medial collagen constituent"""
    return pres_prefactor(x, params) * sigma_collagen_me(x, params)


def pressure_collagen_ad(x, params):
    """Pressure from adventitial collagen constituent"""
    return pres_prefactor(x, params) * sigma_collagen_ad(x, params)


# Force balance equation (solves systolic stretch)
def force_balance_equation(lambda_sys_guess, mE_M, mC_M, mC_A, mM, params):
    """Residual of the force balance equation for root-finding via scipy.fsolve.

    Parameters:
    lambda_sys_guess : sequence of float
        Single-element guess for the systolic stretch (passed by ``fsolve``).
    mE_M, mC_M, mC_A, mM : float
        Current relative mass densities of medial elastin, medial collagen,
        adventitial collagen, and smooth muscle, respectively.
    params : ArterialParameters
        Model parameters.

    Returns:
    float
        Residual of the force balance equation (calculated pressure - target pressure).
    """
    lambda_sys = lambda_sys_guess[0]

    # Individual stresses
    stress_elastin = sigma_elastin(lambda_sys, params)
    stress_collagen_me = sigma_collagen_me(lambda_sys, params)
    stress_collagen_ad = sigma_collagen_ad(lambda_sys, params)
    stress_muscle = sigma_muscle_t(lambda_sys, params)
    mass_density_mult = (
        params.thickness_me * (mE_M * stress_elastin)
        + params.thickness_me * (mC_M * stress_collagen_me)
        + params.thickness_ad * (mC_A * stress_collagen_ad)
        + params.thickness_me * (mM * stress_muscle)
    )
    calculated_pressure = (
        1 / (params.radius_tzero * lambda_sys**2 * params.lambda_z)
    ) * mass_density_mult

    # Return residual from target pressure
    return calculated_pressure - params.pressure_sys


# Medial degeneration by immune cell infiltration
def calculate_immune_cell_level(t, params):
    """
    Immune cell level as a function of time.

    Level is constant at i_0 until time t_i0, after which it increases according
    to a Michaelis-Menten function.
    """
    if t <= params.t_i0:
        return params.i_0
    else:
        return params.i_0 + params.i_max * (
            (t - params.t_i0) / (params.k_i + (t - params.t_i0))
        )


def d_medial_elastin_dt(elastases, elastin_me, params):
    """
    Medial elastin ODE: Degrades due to elastases
    """
    return -params.r_e * elastases * elastin_me


def d_medial_collagen_dt(collagenases, collagen_me, params):
    """
    Medial collagen ODE: Degrades due to collagenases
    """
    return -params.r_cm * collagenases * collagen_me


def d_collagenases_dt(immune_cells, collagenases, params):
    """
    Collagenases ODE: Secretion by immune cells minus decay.
    """
    return params.r_pc1 * immune_cells - params.r_pc2 * collagenases


def d_elastases_dt(immune_cells, elastases, params):
    """
    Elastases ODE: Secretion by immune cells minus decay.
    """
    return params.r_pc1 * immune_cells - params.r_pc2 * elastases


# Adventitial growth and remodeling
def f_lambda_fibroblast(lambda_c_max, lambda_att_max):
    """
    Fibroblast mechanotransduction signal f(lambda).

    Fractional overstretch of collagen beyond its maximum attachment stretch,
    stimulates latent TGF-beta production.
    """
    if lambda_c_max <= lambda_att_max:
        return 0
    return (lambda_c_max - lambda_att_max) / lambda_att_max


def d_fibroblast_dt(tgf_beta, fibroblast, params):
    """
    Fibroblast ODE: TGF-beta-stimulated proliferation minus cell death.

    Fibroblasts are the cells responsible for producing collagen and elastin in the ECM.
    """
    return (
        params.r_f1 + params.r_f2 * tgf_beta
    ) * fibroblast - params.r_f3 * fibroblast


def d_procollagen_dt(tgf_beta, fibroblast, procollagen, params):
    """
    Procollagen ODE: fibroblast secretion minus decay.

    Procollagen is the precursor to collagen, activated by collagenase and
    later matured into collagen.
    """
    return (
        params.r_p1 + params.r_p2 * tgf_beta
    ) * fibroblast - params.r_p3 * procollagen


def d_collagen_dt(procollagen, collagenase, collagen, params):
    """Adventitial collagen ODE: maturation minus collagenase degradation."""
    return params.r_c1 * procollagen - params.r_c2 * collagenase * collagen


def d_zymogen_dt(tgf_beta, fibroblast, zymogen, params):
    """
    Zymogen ODE: fibroblast secretion minus decay.

    Zymogen is the inactive precursor to collagenase, which is activated by TIMPs.
    """
    return (
        params.r_z1 / (1 + params.r_z2 * tgf_beta)
    ) * fibroblast - params.r_z3 * zymogen


def d_collagenase_dt(collagenase, zymogen, timp, params):
    """
    Collagenase ODE: activated by zymogen and inhibited by TIMPs.
    """
    return params.r_ca1 * zymogen - (params.r_ca2 + params.r_ca3 * timp) * collagenase


def d_timp_dt(tgf_beta, fibroblast, collagenase, timp, params):
    """
    TIMP ODE: secretion by fibroblasts and decay/inhibition by collagenase.
    TIMP is the tissue inhibitor of metalloproteinases, which inhibits collagenase.
    """
    return (params.r_i1 + params.r_i2 * tgf_beta) * fibroblast - (
        params.r_i3 + params.r_i4 * collagenase
    ) * timp


def d_latent_tgf_beta_dt(
    tgf_beta,
    latent_tgf_beta,
    fibroblast,
    collagen,
    lambda_c_max,
    lambda_att_max,
    params,
):
    """Latent TGF-beta ODE (Eq. 16): strain-driven production minus turnover.

    Production is governed by the fibroblast mechanotransduction signal and
    collagen level, and scaled by the genotype-dependent baseline
    ``params.tgf_beta_level``.
    """
    term1 = params.r_betal1 * tgf_beta + params.r_betal2 * f_lambda_fibroblast(
        lambda_c_max, lambda_att_max
    )
    term2 = 1 + params.r_betal3 * collagen
    term3 = (
        params.r_betal4
        + params.r_betal5
        * f_lambda_fibroblast(lambda_c_max, lambda_att_max)
        * fibroblast
    )
    # print(f"Bl1 term: {params.r_betal1 * tgf_beta}, Bl2 term: {params.r_betal2 * f_lambda_fibroblast(lambda_c_max, lambda_att_max)}, Bl3 term: {1 + params.r_betal3 * collagen},  Bl4 term: {params.r_betal4}, Bl5 term: {params.r_betal5 * f_lambda_fibroblast(lambda_c_max, lambda_att_max) * fibroblast}")
    return (
        term1 / term2
    ) * fibroblast * params.tgf_beta_level - term3 * latent_tgf_beta


def d_active_tgf_beta_dt(
    tgf_beta, latent_tgf_beta, fibroblast, lambda_c_max, lambda_att_max, params
):
    """
    Active TGF-beta ODE: activation from latent pool minus decay.

    Active TGF-beta is the active form of TGF-beta,
    which stimulates fibroblast proliferation and collagen synthesis.
    """
    return (
        params.r_beta1
        + params.r_beta2
        * f_lambda_fibroblast(lambda_c_max, lambda_att_max)
        * fibroblast
    ) * latent_tgf_beta - params.r_beta3 * tgf_beta


def d_muscle_cells_dt(x, muscle_cells, elastin_me, immune_cells, params):
    """Smooth muscle cell ODE.

    vSMCs proliferate with deviation from homeostatic stretch and are lost as immune
    cells rise; an elastin-coupling term is available but weighted by
    ``beta2_smc`` (currently 0).
    """
    epsilon_stretch = max(0, (x - params.lambda_sys) / params.lambda_sys)
    epsilon_elastin = (elastin_me - params.init_elastin_me) / params.init_elastin_me
    epsilon_immune = (params.i_0 - immune_cells) / 1.0
    return muscle_cells * (
        params.beta1_smc * epsilon_stretch
        + params.beta2_smc * epsilon_elastin
        + params.beta3_smc * epsilon_immune
    )


def calculate_max_attachment_stretch(lambda_c_max_history, dt, t_idx, params):
    """Maximum attachment stretch as a moving average over the remodelling window.

    Averages the historical maximum collagen stretch over the remodelling time
    window (N steps). Ghost values are filled with the initial maximum collagen
    stretch to take into account the initial state of the tissue.
    """
    N = int(params.remodel_time / dt)
    baseline = lambda_c_max_history[0]

    if t_idx < N:
        missing = N - (t_idx + 1)
        lambda_att_max = (
            missing * baseline + np.sum(lambda_c_max_history[: t_idx + 1])
        ) / N
    else:
        window = lambda_c_max_history[t_idx - N + 1 : t_idx + 1]
        lambda_att_max = np.mean(window)

    return lambda_att_max


def calculate_min_attachment_stretch(lambda_att_max, params):
    """
    Minimum attachment stretch: max minus the distribution width.

    Width w(t) is taken as constant, equal to ``params.width_att_dist``.
    """
    lambda_att_min = lambda_att_max - params.width_att_dist
    return lambda_att_min


def calculate_mode_attachment_stretch(lambda_att_min, lambda_att_max, params):
    """
    Mode attachment stretch: min plus skew * (max - min)

    The skew s(t) is taken as constant, equal to ``params.skew_att_dist``.
    """
    lambda_att_mode = lambda_att_min + params.skew_att_dist * (
        lambda_att_max - lambda_att_min
    )
    return lambda_att_mode


def calculate_min_attachment_stretch_me(lambda_att_max_me, params):
    """Medial minimum attachment stretch (max minus medial distribution width)."""
    return lambda_att_max_me - params.width_att_dist_me


def calculate_mode_attachment_stretch_me(lambda_att_min_me, lambda_att_max_me, params):
    """Medial mode attachment stretch (min plus medial skew * (max - min))."""
    return lambda_att_min_me + params.skew_att_dist_me * (
        lambda_att_max_me - lambda_att_min_me
    )


# Collagen recruitment stretch remodelling equations
def alpha_rate(fibroblast, collagen, collagenase, params):
    """
    Fibroblast remodelling rate alpha

    Balances fibroblast activity against collagen availability, with a
    synthesis/degradation term that scales with the square root of collagen
    and collagenase levels.
    """
    alpha = (
        params.alpha_init * (fibroblast / collagen) * np.sqrt(collagen * collagenase)
    )
    return alpha


def d_collagen_min_recruitment_stretch_ad_dt(alpha, lambda_c_max, lambda_att_max):
    """
    Rate of change of the adventitial minimum recruitment stretch
    """
    return alpha * (lambda_c_max - lambda_att_max) / lambda_att_max


def d_collagen_max_recruitment_stretch_ad_dt(alpha, lambda_c_min, lambda_att_min):
    """
    Rate of change of the adventitial maximum recruitment stretch
    """
    return alpha * (lambda_c_min - lambda_att_min) / lambda_att_min


def d_collagen_mode_recruitment_stretch_ad_dt(alpha, lambda_c_mode, lambda_att_mode):
    """
    Rate of change of the mode recruitment stretch for collagen in the adventitia.
    """
    return alpha * (lambda_c_mode - lambda_att_mode) / lambda_att_mode


# Add TGF-beta1 protein level function
def get_latent_tgf_beta_level(params, genotype=None):
    """Return the genotype-specific baseline latent TGF-beta level.

    Looks up ``genotype`` in ``params.tgf_beta_levels``, returns
    1.0 if none or Nan. Used to scale latent TGF-beta production.

    Parameters
    params : ArterialParameters
        Provides the genotype-to-level lookup table.
    genotype : str, optional
        Patient genotype, e.g. "TT", "TC", "CC" (case-insensitive).
    """
    levels = params.tgf_beta_levels
    if genotype is None:
        return 1.0
    return levels.get(genotype.upper(), 1.0)
