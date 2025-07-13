import numpy as np

class ArterialParameters:
    def __init__(self):
        # Geometric and pressure
        self.c_diam_tzero_mm = 2.9
        self.c_radius_tzero = (self.c_diam_tzero_mm * 1e-3) / (2 * 1.3)
        self.c_thickness_tzero = self.c_radius_tzero / 5
        self.c_pressure_sys = 16000

        # Stretches
        self.c_lambda_z = 1.3
        self.c_lambda_sys = 1.3 
        self.c_lambda_elastin = 1.3 
        self.c_lambda_muscle = 1.15
        self.c_rec_muscle = self.c_lambda_elastin / self.c_lambda_muscle
        self.c_musc_mean = 1.1
        self.c_musc_min = 0.4
        self.c_vasodil_conc = 0.68
        self.c_ge_muscle = (self.c_lambda_muscle**2 - 1.0) / 2.0

        # Collagen ratio
        self.c_collagen_ratio_ad_me = 8

        # Media attachment stretches
        self.c_att_min_me = 1.00001
        self.c_att_mod_me = 1.01
        self.c_att_max_me = 1.07

        # Media recruitment stretches
        self.c_rec_max_me = self.c_lambda_elastin / self.c_att_min_me
        self.c_rec_min_me = self.c_lambda_elastin / self.c_att_max_me
        self.c_rec_mod_me = self.c_lambda_elastin / self.c_att_mod_me
        self.v_a_me = self.c_rec_min_me
        self.v_b_me = self.c_rec_max_me
        self.v_c_me = self.c_rec_mod_me

        # Adventitia attachment
        self.c_att_min_ad = 0.8
        self.c_att_mod_ad = 0.9
        self.c_att_max_ad = 0.99999

        # Adventitia recruitment
        self.c_rec_max_ad = self.c_lambda_elastin / self.c_att_min_ad
        self.c_rec_min_ad = self.c_lambda_elastin / self.c_att_max_ad
        self.c_rec_mod_ad = self.c_lambda_elastin / self.c_att_mod_ad
        self.v_a_ad = self.c_rec_min_ad
        self.v_b_ad = self.c_rec_max_ad
        self.v_c_ad = self.c_rec_mod_ad

        # Load-borne proportions
        self.c_load_borne_elastin = 0.50
        self.c_load_borne_muscle_p = 0.20
        self.c_load_borne_muscle_a = 0.20
        self.c_load_borne_collagen = 1.0 - (
            self.c_load_borne_elastin
            + self.c_load_borne_muscle_p
            + self.c_load_borne_muscle_a
        )

        # Common factor for tension balance
        self.c_common_factor = (
            self.c_pressure_sys
            * self.c_radius_tzero
            * self.c_lambda_elastin**2
            * self.c_lambda_z
            / self.c_thickness_tzero
        )

        # Material parameters
        self.c_k_elastin = (
            self.c_load_borne_elastin
            * self.c_common_factor
            / (self.c_lambda_elastin**2 * (1 - (1 / (self.c_lambda_z**2 * self.c_lambda_elastin**4))))
        ) 

        collagen_denominator = (
            2
            * self.c_lambda_elastin
            / ((self.v_b_me - self.v_a_me) * (self.v_c_me - self.v_a_me))
            * ((self.v_a_me + self.c_lambda_elastin) * np.log(self.c_lambda_elastin / self.v_a_me)
               + 2 * (self.v_a_me - self.c_lambda_elastin))
        )

        self.c_k_collagen = self.c_load_borne_collagen * self.c_common_factor / collagen_denominator

        # Media collagen Cauchy stress
        self.v_gamma_me = self.c_k_collagen / ((self.v_b_me - self.v_a_me) * (self.v_c_me - self.v_a_me))
        self.v_delta_me = self.c_k_collagen / ((self.v_b_me - self.v_a_me) * (self.v_b_me - self.v_c_me))

        # Adventitia collagen Cauchy stress
        self.v_gamma_ad = self.c_k_collagen * self.c_collagen_ratio_ad_me / ((self.v_b_ad - self.v_a_ad) * (self.v_c_ad - self.v_a_ad))
        self.v_delta_ad = self.c_k_collagen * self.c_collagen_ratio_ad_me / ((self.v_b_ad - self.v_a_ad) * (self.v_b_ad - self.v_c_ad))


        muscle_a_denominator = (
            self.c_vasodil_conc
            * self.c_lambda_muscle
            * (1 - ((self.c_musc_mean - self.c_lambda_muscle) / (self.c_musc_mean - self.c_musc_min)) ** 2)
        )

        self.c_k_muscle_p = (
            self.c_load_borne_muscle_p
            * self.c_common_factor
            / (self.c_lambda_muscle**2 * (1 - 1 / (self.c_lambda_z**2 * self.c_lambda_muscle**4)))
        )

        self.c_k_muscle_a = self.c_load_borne_muscle_a * self.c_common_factor / muscle_a_denominator

        # Immune cell related rates
        self.r_e = 1.0       # Elastin degradation rate by immune cell proteases (years^-1)
        self.r_cm = 1.0      # Medial collagen degradation rate by immune cell proteases (years^-1)
        self.t_i0 = 40       # Time for onset of immune cell infiltration (years)
        self.i_0 = 0         # Initial level of immune cells in the arterial wall
        self.i_max = 1.0     # Maximum level of immune cells in the arterial wall
        self.k_i = 1.25      # Time for immune cell levels to reach half of maximum level (years)
        self.r_pc1 = 1.0     # Collagenase secretion rate by immune cells (years^-1)
        self.r_pc2 = 1.0     # Baseline immune cell collagenase degradation rate (years^-1)
        self.r_pe1 = 1.0     # Elastase secretion rate by immune cells (years^-1)
        self.r_pe2 = 1.0     # Baseline immune cell elastase degradation rate (years^-1)

        # Fibroblast rates
        self.r_f1 = 1.0      # Baseline fibroblast migration and proliferation rate (years^-1)
        self.r_f2 = 0.5      # Fibroblast population dynamics sensitivity to TGF-Beta (years^-1)
        self.r_f3 = 1.0      # Fibroblast cell death rate (years^-1)

        # Procollagen rates
        self.r_p1 = 1.0      # Baseline procollagen secretion rate by fibroblasts (years^-1)
        self.r_p2 = 0.5      # Fibroblast procollagen secretion sensitivity to TGF-Beta (years^-1)
        self.r_p3 = 1.0      # Combined baseline procollagen degradation and modification rate (years^-1)

        # Adventitial collagen rates
        self.r_c1 = 0.5      # Baseline adventitial collagen maturation rat (years^-1)
        self.r_c2 = 0.5      # Adventitial collagen degredation rate (years^-1)

        # Zymogen rates
        self.r_z1 = 1.0      # Baseline zymogen secretion rate by fibroblasts (years^-1)
        self.r_z2 = 0.5      # Fibroblast zymogen secretion sensitivity to TGF-Beta (years^-1)
        self.r_z3 = 1.0      # Combined baseline zymogen degradation and modification rates (years^-1)

        # Collagenase rates
        self.r_ca1 = 0.5      # Baseline active collagenase maturation rate (years^-1)
        self.r_ca2 = 0.25    # Collagenase inactivation rate by TIMPs (years^-1)
        self.r_ca3 = 0.25    # Inhibitor–collagenase complex formation rate (years^-1)

        # Inhibitor rates
        self.r_i1 = 1.0       # Baseline inhibitor secretion rate (years^-1)
        self.r_i2 = 0.5       # Fibroblast collagenase inhibitor secretion sensitivity to TGF-β (years^-1)
        self.r_i3 = 0.75      # Baseline collagenase inhibitor degradation rate (years^-1)
        self.r_i4 = 0.25      # Inhibitor–collagenase complex formation rate (duplicate of r_ca3) (years^-1)

        # TGF-Beta rates
        self.r_betal1 = 0.1          # Fibroblast latent TGF-Beta secretion sensitivity to active TGF-Beta
        self.r_betal2 = 5.0          # Fibroblast latent TGF-Beta secretion sensitivity to deviations from mechanical homeostasis (Parameter study with [0.1, 1.0, 5.0, 10.0])
        self.r_betal3 = 1.0          # Fibroblast latent TGF-Beta secretion sensitivity to collagen levels
        self.r_betal4 = 1.0          # Combined baseline latent TGF-Bet degradation/modification rate (years^-1)
        self.r_betal5 = 1.0          # Latent TGF-Beta modification rate by integrin/ECM/Stretch–dependent mechanism (years^-1
        self.r_beta1 = 0.5           # Baseline latent TGF-Beta activation rate
        self.r_beta2 = 1.0           # Modification rate by integrin/ECM/Stretch–dependent mechanism (same as r_betaL5)
        self.r_beta3 = 1.0           # Baseline active TGF-Beta degradation rate (years^-1)

        # Mechanical model parameters
        self.remodel_time = 10       # Averaging time period for attachment stretch remodelling 
        self.t_sim = 90              # Simulation time in years
        self.width_att_dist = 0.1    # Width of the attachment stretch distribution (assumed constant)
        self.skew_att_dist = 0.5     # Skew of the attachment stretch distribution (assumed constant)

        # Set initial values for variables 
        self.alpha_init = 1.15
        self.init_fibroblast = 1.0
        self.init_collagen_ad = 1.0
        self.init_collagen_me = 1.0
        self.init_elastin_ad = 1.0
        self.init_elastin_me = 1.0
        self.init_procollagen = 1.0
        self.init_collagenase = 1.0
        self.init_zymogen = 1.0
        self.init_timp = 1.0
        self.init_collagenases = 0.0
        self.init_elastases = 0.0
        self.init_latent_tgf_beta = 1e-3
        self.init_active_tgf_beta = 1.0


    def to_dict(self):
        return self.__dict__
