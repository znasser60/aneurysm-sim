import numpy as np
import matplotlib.pyplot as plt

from config import ArterialParameters

def main(): 
    """
    Main function to run the arterial parameter extraction and display results.
    """
    params = ArterialParameters()
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
    
    # Plotting the results
    n_zoom = 120

    plt.figure(figsize=(12, 8))
    plt.plot(sv_stretch_var[:n_zoom], sv_pressure_var[:n_zoom]/1e3, linewidth=2, label='Total')
    plt.plot(sv_stretch_var[32:n_zoom], sv_pressure_var_elastin[32:n_zoom]/1e3, '--', linewidth=2, label='Elastin')
    plt.plot(sv_stretch_var[64:n_zoom], sv_pressure_var_collagen[64:n_zoom]/1e3, '-.', linewidth=2, label='Collagen')
    plt.plot(sv_stretch_var[64:n_zoom], sv_pressure_var_collagen_me[64:n_zoom]/1e3, '--', linewidth=2, label='Collagen Media')
    plt.plot(sv_stretch_var[64:n_zoom], sv_pressure_var_collagen_ad[64:n_zoom]/1e3, '--', linewidth=2, label='Collagen Adventitia')
    plt.plot(sv_stretch_var[44:n_zoom], sv_pressure_var_muscle_p[44:n_zoom]/1e3, '--', linewidth=2, label='Muscle Passive')
    plt.plot(sv_stretch_var[:n_zoom], sv_pressure_var_muscle_a[:n_zoom]/1e3, '--', linewidth=2, label='Muscle Active')
    plt.title('Pressure vs Stretch')
    plt.xlabel('Stretch')
    plt.ylabel('Pressure (kPa)')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.ylim([0, 60])
    plt.show()

    # Plot 2: Pressure vs Diameter 
    n_zoom_diam = 125
    sv_diam_var = 2 * params.c_radius_tzero * 1e3 * sv_stretch_var
    sv_pressure_var_smc = sv_pressure_var_muscle_a + sv_pressure_var_muscle_p
    sv_pressure_var_coll = sv_pressure_var_collagen_me + sv_pressure_var_collagen_ad
    
    plt.figure(figsize=(12, 8))
    plt.plot(sv_diam_var[:n_zoom_diam], sv_pressure_var[:n_zoom_diam]/1e3, 'k-', linewidth=6, label='Total')
    plt.plot(sv_diam_var[32:n_zoom_diam], sv_pressure_var_elastin[32:n_zoom_diam]/1e3, 'k--', linewidth=3, label='E')
    plt.plot(sv_diam_var[66:n_zoom_diam], sv_pressure_var_coll[66:n_zoom_diam]/1e3, 'k-.', linewidth=3, label='C')
    plt.plot(sv_diam_var[44:n_zoom_diam], sv_pressure_var_muscle_p[44:n_zoom_diam]/1e3, 'k+', linewidth=2, markersize=8, label='VSMCp')
    plt.plot(sv_diam_var[:n_zoom_diam], sv_pressure_var_muscle_a[:n_zoom_diam]/1e3, 'k.', linewidth=2, markersize=8, label='VSMCa')
    
    plt.axvline(x=2.9, color='red', linestyle='--', linewidth=2)
    plt.axhline(y=16, color='red', linestyle='--', linewidth=2)
    plt.plot(2.9, 16, 'ro', linewidth=3, markersize=8)
    
    plt.title('Pressure vs Diameter')
    plt.xlabel('Diameter (mm)')
    plt.ylabel('Pressure (kPa)')
    plt.legend(loc='upper left')
    plt.ylim([0, 60])
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()
