import numpy as np
from scipy.optimize import fsolve
 
from aneurysm_sim.model import functions
from aneurysm_sim.config.parameters import ArterialParameters
 
 
def simulate_arterial_stress_and_pressure(params):
    """Computes individual component stress and pressure curves over a circumferential stretch sweep.


    Parameters
    params : ArterialParameters
        Model parameters for a single patient/genotype.
 
    Returns
    dict
        Stretch and diameter axes plus constituent stress and pressure arrays.
    """
    n = 235
    stretch = 0.55 + 0.01 * np.arange(n)  # circumferential stretch sweep
    diam = 2 * params.radius_tzero * stretch
 
    # Stress curves (scalar constitutive laws → looped).
    stress = {k: np.zeros(n) for k in (
        "elastin", "collagen_me", "collagen_ad", "collagen",
        "muscle_a", "muscle_p", "muscle_t", "total",
    )}
    pressure = {k: np.zeros(n) for k in (
        "total", "elastin", "collagen", "muscle", "muscle_a", "muscle_p",
        "collagen_me", "collagen_ad",
    )}
 
    for i in range(n):
        x = stretch[i]
 
        stress["elastin"][i] = functions.sigma_elastin(x, params)
        stress["collagen"][i] = functions.sigma_collagen(x, params)
        stress["collagen_me"][i] = functions.sigma_collagen_me(x, params)
        stress["collagen_ad"][i] = functions.sigma_collagen_ad(x, params)
        stress["muscle_a"][i] = functions.sigma_muscle_a(x, params)
        stress["muscle_p"][i] = functions.sigma_muscle_p(x, params)
        stress["muscle_t"][i] = functions.sigma_muscle_t(x, params)
        stress["total"][i] = max(
            stress["elastin"][i] + stress["collagen"][i] + stress["muscle_t"][i], 0
        )
 
        pressure["elastin"][i] = max(functions.pressure_elastin(x, params), 0)
        pressure["collagen"][i] = functions.pressure_collagen(x, params)
        pressure["collagen_me"][i] = functions.pressure_collagen_me(x, params)
        pressure["collagen_ad"][i] = functions.pressure_collagen_ad(x, params)
        pressure["muscle_p"][i] = max(functions.pressure_muscle_p(x, params), 0)
        pressure["muscle_a"][i] = max(functions.pressure_muscle_a(x, params), 0)
        pressure["muscle"][i] = pressure["muscle_a"][i] + pressure["muscle_p"][i]
        pressure["total"][i] = (
            pressure["elastin"][i]
            + pressure["collagen_me"][i]
            + pressure["collagen_ad"][i]
            + pressure["muscle"][i]
        )
 
    out = {"sv_diam_var": diam, "sv_stretch_var": stretch}
    out.update({f"sv_stress_var_{k}": v for k, v in stress.items()})
    out["sv_pressure_var"] = pressure.pop("total")
    out.update({f"sv_pressure_var_{k}": v for k, v in pressure.items()})
    return out
 
 
def _apply_collagen_shape_params(params, layer, rec_min, rec_max, rec_mode, stiffness_scale):
    """Update collagen distribution parameters based on the current recruitment-stretch evolution.
 
    Parameters
    params : ArterialParameters
        Parameter object mutated in place.
    layer : {"me", "ad"}
        Layer suffix to update.
    rec_min, rec_max, rec_mode : float
        Current minimum/maximum/mode recruitment stretches.
    stiffness_scale : float
        Collagen stiffness multiplier (1.0 for media, the adventitia:media
        collagen ratio (8.0) for adventitia).
    """
    a, b, c = rec_min, rec_max, rec_mode
    setattr(params, f"a_{layer}", a)
    setattr(params, f"b_{layer}", b)
    setattr(params, f"c_{layer}", c)
    k = params.k_collagen * stiffness_scale
    setattr(params, f"gamma_{layer}", k / ((b - a) * (c - a)))
    setattr(params, f"delta_{layer}", k / ((b - a) * (b - c)))
 
 
def simulate_aneurysm(params, treatment=False, dt=0.01):
    """Integrate the aneurysm growth and remodelling model over time.
 
    Holds the artery at a healthy homeostasis until the onset of the immune cell infiltration 
    then integrates the ODEs (see ``functions.py``) for the medial/adventitial remodelling and the 
    force balance to determine the systolic stretch.  
 
    Parameters
    params : ArterialParameters
        Model parameters for the patient/genotype being simulated.
    treatment : bool, default False
        If True apply the single active-TGF-beta pulse at treatment time ``params.t_treat``.
    dt : float, default 0.01
        Integration time step in years (found via a convergence study).
 
    Returns
    -------
    dict
        Time, mass-density trajectories, collagen attachment/recruitment-stretch histories, 
        the systolic-stretch history, and the final systolic stretch.
    """
    steps = int(params.t_sim / dt) + 1
    time = np.linspace(0, params.t_sim, steps)
 
    def _traj():
        # Helper to create a zero trajectory array of the correct length
        return np.zeros(steps)

    # INITIALIZE TRAJECTORIES 
    # Medial degradation variables
    elastin_me, collagen_me = _traj(), _traj()
    elastases, collagenases, immune_cells = _traj(), _traj(), _traj()
    muscle_cells = _traj()
    # Adventitial remodelling variables
    collagen_ad, fibroblast, procollagen = _traj(), _traj(), _traj()
    collagenase, zymogen, timp = _traj(), _traj(), _traj()
    latent_tgf_beta, active_tgf_beta = _traj(), _traj()
    # Systolic stretch, diameter, adventitial collagen stretch variables
    diameter, lambda_sys_array = _traj(), _traj()
    lambda_c_max, lambda_c_min, lambda_c_mode = _traj(), _traj(), _traj()
    lambda_att_max, lambda_att_min, lambda_att_mode = _traj(), _traj(), _traj()
    lambda_rec_max, lambda_rec_min, lambda_rec_mode = _traj(), _traj(), _traj()
    # Medial collagen stretch variables
    lambda_c_max_me, lambda_c_min_me, lambda_c_mode_me = _traj(), _traj(), _traj()
    lambda_att_max_me, lambda_att_min_me, lambda_att_mode_me = _traj(), _traj(), _traj()
    lambda_rec_min_me, lambda_rec_max_me, lambda_rec_mode_me = _traj(), _traj(), _traj()
 
    # Homeostatic values
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
 
    lambda_sys_array[0] = params.lambda_sys
    diameter[0] = 2 * params.radius_tzero * params.lambda_sys
 
    # Adventitial collagen attachment/recruitment distributions
    lambda_att_max[0] = params.att_max_ad
    lambda_att_min[0] = params.att_min_ad
    lambda_att_mode[0] = params.att_mod_ad
    lambda_rec_max[0] = params.rec_max_ad
    lambda_rec_min[0] = params.rec_min_ad
    lambda_rec_mode[0] = params.rec_mod_ad
    lambda_c_max[0] = params.lambda_sys / lambda_rec_min[0]
    lambda_c_min[0] = params.lambda_sys / lambda_rec_max[0]
    lambda_c_mode[0] = params.lambda_sys / lambda_rec_mode[0]
 
    # Medial collagen attachment/recruitment distributions
    lambda_att_max_me[0] = params.att_max_me
    lambda_att_min_me[0] = params.att_min_me
    lambda_att_mode_me[0] = params.att_mod_me
    lambda_rec_min_me[0] = params.rec_min_me
    lambda_rec_max_me[0] = params.rec_max_me
    lambda_rec_mode_me[0] = params.rec_mod_me
    lambda_c_max_me[0] = params.lambda_sys / lambda_rec_min_me[0]
    lambda_c_min_me[0] = params.lambda_sys / lambda_rec_max_me[0]
    lambda_c_mode_me[0] = params.lambda_sys / lambda_rec_mode_me[0]
 
    # Hold the artery at homeostasis until the immune cell infiltration onset
    onset_idx = int(np.count_nonzero(time < params.t_i0))
    held = [
        elastin_me, collagen_me, elastases, collagenases, immune_cells,
        muscle_cells, collagen_ad, fibroblast, procollagen, collagenase,
        zymogen, timp, latent_tgf_beta, active_tgf_beta,
        diameter, lambda_sys_array,
        lambda_c_max, lambda_c_min, lambda_c_mode,
        lambda_att_max, lambda_att_min, lambda_att_mode,
        lambda_rec_max, lambda_rec_min, lambda_rec_mode,
        lambda_c_max_me, lambda_c_min_me, lambda_c_mode_me,
        lambda_att_max_me, lambda_att_min_me, lambda_att_mode_me,
        lambda_rec_min_me, lambda_rec_max_me, lambda_rec_mode_me,
    ]
    for arr in held:
        arr[:onset_idx] = arr[0]
    lambd_c_max_history = [lambda_c_max[0]] * onset_idx
    lambd_c_max_me_history = [lambda_c_max_me[0]] * onset_idx
 
    # Integrate the ODEs forward in time from the onset of immune cell infiltration
    for i in range(onset_idx, steps):
        t = time[i]
 
        # Biological ODEs (fwd Euler integration)
        immune_cells[i] = functions.calculate_immune_cell_level(t, params)
        collagen_me[i] = collagen_me[i - 1] + dt * functions.d_medial_collagen_dt(
            collagenases[i - 1], collagen_me[i - 1], params
        )
        collagen_ad[i] = collagen_ad[i - 1] + dt * functions.d_collagen_dt(
            procollagen[i - 1], collagenase[i - 1], collagen_ad[i - 1], params
        )
        elastin_me[i] = elastin_me[i - 1] + dt * functions.d_medial_elastin_dt(
            elastases[i - 1], elastin_me[i - 1], params
        )
        collagenases[i] = collagenases[i - 1] + dt * functions.d_collagenases_dt(
            immune_cells[i - 1], collagenases[i - 1], params
        )
        elastases[i] = elastases[i - 1] + dt * functions.d_elastases_dt(
            immune_cells[i - 1], elastases[i - 1], params
        )
        muscle_cells[i] = muscle_cells[i - 1] + dt * functions.d_muscle_cells_dt(
            lambda_sys_array[i - 1], muscle_cells[i - 1], elastin_me[i - 1],
            immune_cells[i - 1], params,
        )
        fibroblast[i] = fibroblast[i - 1] + dt * functions.d_fibroblast_dt(
            active_tgf_beta[i - 1], fibroblast[i - 1], params
        )
        procollagen[i] = procollagen[i - 1] + dt * functions.d_procollagen_dt(
            active_tgf_beta[i - 1], fibroblast[i - 1], procollagen[i - 1], params
        )
        collagenase[i] = collagenase[i - 1] + dt * functions.d_collagenase_dt(
            collagenase[i - 1], zymogen[i - 1], timp[i - 1], params
        )
        zymogen[i] = zymogen[i - 1] + dt * functions.d_zymogen_dt(
            active_tgf_beta[i - 1], fibroblast[i - 1], zymogen[i - 1], params
        )
        timp[i] = timp[i - 1] + dt * functions.d_timp_dt(
            active_tgf_beta[i - 1], fibroblast[i - 1], collagenase[i - 1],
            timp[i - 1], params,
        )
        latent_tgf_beta[i] = latent_tgf_beta[i - 1] + dt * functions.d_latent_tgf_beta_dt(
            active_tgf_beta[i - 1], latent_tgf_beta[i - 1], fibroblast[i - 1],
            collagen_ad[i - 1], lambda_c_max[i - 1], lambda_att_max[i - 1], params,
        )
        active_tgf_beta[i] = active_tgf_beta[i - 1] + dt * functions.d_active_tgf_beta_dt(
            active_tgf_beta[i - 1], latent_tgf_beta[i - 1], fibroblast[i - 1],
            lambda_c_max[i - 1], lambda_att_max[i - 1], params,
        )
 
        # Adventitial recruitment-stretch evolution 
        alpha = functions.alpha_rate(fibroblast[i], collagen_ad[i], collagenase[i], params)
        lambda_rec_min[i] = lambda_rec_min[i - 1] + dt * functions.d_collagen_min_recruitment_stretch_ad_dt(
            alpha, lambda_c_max[i - 1], lambda_att_max[i - 1]
        )
        lambda_rec_max[i] = lambda_rec_max[i - 1] + dt * functions.d_collagen_max_recruitment_stretch_ad_dt(
            alpha, lambda_c_min[i - 1], lambda_att_min[i - 1]
        )
        lambda_rec_mode[i] = lambda_rec_mode[i - 1] + dt * functions.d_collagen_mode_recruitment_stretch_ad_dt(
            alpha, lambda_c_mode[i - 1], lambda_att_mode[i - 1]
        )
        _apply_collagen_shape_params(
            params, "ad", lambda_rec_min[i], lambda_rec_max[i], lambda_rec_mode[i],
            params.collagen_ratio_ad_me,
        )
 
        # Medial recruitment-stretch evolution 
        lambda_rec_min_me[i] = lambda_rec_min_me[i - 1] + dt * functions.d_collagen_min_recruitment_stretch_ad_dt(
            params.alpha_init, lambda_c_max_me[i - 1], lambda_att_max_me[i - 1]
        )
        lambda_rec_max_me[i] = lambda_rec_max_me[i - 1] + dt * functions.d_collagen_max_recruitment_stretch_ad_dt(
            params.alpha_init, lambda_c_min_me[i - 1], lambda_att_min_me[i - 1]
        )
        lambda_rec_mode_me[i] = lambda_rec_mode_me[i - 1] + dt * functions.d_collagen_mode_recruitment_stretch_ad_dt(
            params.alpha_init, lambda_c_mode_me[i - 1], lambda_att_mode_me[i - 1]
        )
        _apply_collagen_shape_params(
            params, "me", lambda_rec_min_me[i], lambda_rec_max_me[i], lambda_rec_mode_me[i],
            1.0,
        )
 
        # Solve the force balance eq. for the current systolic stretch
        lambda_sys = fsolve(
            functions.force_balance_equation,
            [lambda_sys_array[i - 1]],
            args=(elastin_me[i], collagen_me[i], collagen_ad[i], muscle_cells[i], params),
        )[0]
        lambda_sys_array[i] = lambda_sys
        diameter[i] = 2 * params.radius_tzero * lambda_sys
 
        # Update adventitial collagen stretches and attachment distribution
        lambda_c_max[i] = lambda_sys / lambda_rec_min[i]
        lambda_c_min[i] = lambda_sys / lambda_rec_max[i]
        lambda_c_mode[i] = lambda_sys / lambda_rec_mode[i]
        lambd_c_max_history.append(lambda_c_max[i])
        lambda_att_max[i] = functions.calculate_max_attachment_stretch(lambd_c_max_history, dt, i, params)
        lambda_att_min[i] = functions.calculate_min_attachment_stretch(lambda_att_max[i], params)
        lambda_att_mode[i] = functions.calculate_mode_attachment_stretch(lambda_att_min[i], lambda_att_max[i], params)
 
        # Update medial collagen stretches and attachment distribution
        lambda_c_max_me[i] = lambda_sys / lambda_rec_min_me[i]
        lambda_c_min_me[i] = lambda_sys / lambda_rec_max_me[i]
        lambda_c_mode_me[i] = lambda_sys / lambda_rec_mode_me[i]
        lambd_c_max_me_history.append(lambda_c_max_me[i])
        lambda_att_max_me[i] = functions.calculate_max_attachment_stretch(lambd_c_max_me_history, dt, i, params)
        lambda_att_min_me[i] = functions.calculate_min_attachment_stretch_me(lambda_att_max_me[i], params)
        lambda_att_mode_me[i] = functions.calculate_mode_attachment_stretch_me(lambda_att_min_me[i], lambda_att_max_me[i], params)
 
        # Optional therapeutic active-TGF-beta pulse
        if treatment and abs(t - params.t_treat) < dt:
            active_tgf_beta[i] += params.tgf_spike_amount
 
    return {
        "diameter": diameter,
        "time": time,
        "elastin_me": elastin_me,
        "collagen_me": collagen_me,
        "elastases": elastases,
        "collagenases": collagenases,
        "immune_cells": immune_cells,
        "fibroblast": fibroblast,
        "muscle_cells": muscle_cells,
        "collagen_ad": collagen_ad,
        "procollagen": procollagen,
        "collagenase": collagenase,
        "zymogen": zymogen,
        "timp": timp,
        "latent_tgf_beta": latent_tgf_beta,
        "active_tgf_beta": active_tgf_beta,
        "lambda_c_max": lambda_c_max,
        "lambda_c_min": lambda_c_min,
        "lambda_c_mode": lambda_c_mode,
        "lambda_att_max": lambda_att_max,
        "lambda_att_min": lambda_att_min,
        "lambda_att_mode": lambda_att_mode,
        "lambda_rec_max": lambda_rec_max,
        "lambda_rec_min": lambda_rec_min,
        "lambda_rec_mode": lambda_rec_mode,
        "lambd_c_max_history": lambd_c_max_history,
        "lambda_sys": lambda_sys_array,
        "final_lambda_sys": lambda_sys_array[-1],
    }
 
 
def simulate_aneurysm_batch_smc(params, treatment=False, n_vals=15):
    """Run a vSMC-fraction sensitivity sweep for one polygenic score.
 
    Samples n values vSMC volume fractions evenly across ±1 SD about the mean
    for the score (clipped to [0, 1])
 
    Parameters
    params : ArterialParameters
        Model parameters for the patient/genotype being simulated
    treatment : bool, default False
        Boolean to apply the therapeutic TGF-beta pulse in each run
    n_vals : int, default 15
        Number of fractions to sample across the ±1 SD range
 
    Returns
    -------
    tuple(list[dict], dict)
        The list of simulation results (per vSMC fraction) and the mean result (for the mean vSMC fraction)
    """
    score = (
        int(np.clip(params.polygenic_score, 0, 4))
        if params.polygenic_score is not None
        else 0
    )
    mean_val = params.smc_mean_fractions[score]
    sd_val = params.smc_sd_fractions[score]
 
    sweep_fractions = np.linspace(
        np.clip(mean_val - sd_val, 0.0, 1.0),
        np.clip(mean_val + sd_val, 0.0, 1.0),
        n_vals,
    )
 
    batch_results = []
    for fraction in sweep_fractions:
        p_temp = ArterialParameters(
            smc_fraction=fraction,
            polygenic_score=params.polygenic_score,
            tgf_beta_level=params.tgf_beta_level,
        )
        try:
            batch_results.append(simulate_aneurysm(p_temp, treatment=treatment))
        except Exception as e:  # keep sweeping if one fraction fails
            print(f"Simulation failed for fraction {fraction:.4f}: {e}")
 
    p_mean = ArterialParameters(
        smc_fraction=mean_val,
        polygenic_score=params.polygenic_score,
        tgf_beta_level=params.tgf_beta_level,
    )
    mean_result = simulate_aneurysm(p_mean, treatment=treatment)
 
    return batch_results, mean_result