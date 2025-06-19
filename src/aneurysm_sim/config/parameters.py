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

    def to_dict(self):
        return self.__dict__
