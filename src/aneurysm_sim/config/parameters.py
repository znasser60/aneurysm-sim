import numpy as np


class ArterialParameters:
    """
    Class to hold general and patient-specific initial parameters for the 1D cerebral artery model.

    A patient is specified by a genotype and/or a polygenic score to initilize
    the simulation. Baseline TGF-beta production level and smooth muscle cell (SMC)
    volume fraction are looked up, and either derived quantity can be overridden
    directly via ```tgf_beta_level``` or ```smc_fraction``` for sensitivity sweeps
    and synthetic-patient studies. The method `refresh_physics` is called at each simulation
    run to calculaate the constituent stiffness constants. 

    Parameters:
    genotype: str, optional
        TGF-B genotype of the patient, one of "TT", "TC", or "CC". If
        none or NaN, the default TGF-beta production level of 1.0 is used.
    polygenic_score: int, optional
        vSMC polygenic score of the patient, an integer between 0 and 4. If
        none or NaN, the default SMC fraction of 0.7262666 (no risk alleles) is used.
    smc_fraction: float, optional
        Override for the smooth muscle cell (SMC) volume fraction in the arterial wall
    tgf_beta_level: float, optional
        Override for the latent TGF-beta production level in the arterial wall


    Notes:
    ________
    The genotype-specific TGF-beta levels and the polygenic-score SMC fraction
    statistics were provided by Dr. Mark Bakker (UMC Utrecht) and are derived
    from UK Biobank pQTL data and bulk RNA-seq cell-type deconvolution,
    respectively.
    """

    def __init__(
        self,
        genotype=None,
        polygenic_score=None,
        smc_fraction=None,
        tgf_beta_level=None,
    ):
        # Geometric and pressure
        self.diam_tzero_mm = 2.9  # Loaded initial diameter in mm
        self.radius_tzero = self.diam_tzero_mm / (
            2 * 1.3
        )  # Loaded initial radius in mm
        self.thickness_ad = 0.104  # Adventitia thickness in mm
        self.thickness_me = 0.216  # Media thickness in mm
        self.thickness_tzero = self.thickness_ad + self.thickness_me
        self.pressure_sys = 16000  # Systolic pressure in Pa (120 mmHg)

        # Stretches
        self.lambda_z = 1.3  # Axial stretch ratio
        self.lambda_sys = 1.3  # Circumferential stretch ratio at systolic pressure
        self.lambda_elastin = 1.3  # Elastin stretch ratio
        self.lambda_muscle = 1.13  # Muscle stretch ratio
        self.lambda_muscle_att = 1.1  # Muscle attachment stretch ratio
        self.rec_muscle = (
            self.lambda_sys / self.lambda_muscle
        )  # Muscle recruitment stretch ratio
        self.musc_mean = 1.1  # Mean muscle stretch ratio
        self.musc_min = 0.4  # Minimum muscle stretch ratio
        self.vasodil_conc = 0.68  # Vasodilator concentration
        self.vasodil_conc_basal = 0.68  # Basal vasodilator concentration
        self.vasodil_conc_shear = 1.36  # Shear-induced vasodilator concentration

        # Collagen adventitia:media ratio
        self.collagen_ratio_ad_me = 8.0

        # Media attachment stretch distribution
        self.att_min_me = 1.00001
        self.att_mod_me = 1.01
        self.att_max_me = 1.07

        # Media recruitment stretch distribution
        self.rec_max_me = self.lambda_elastin / self.att_min_me
        self.rec_min_me = self.lambda_elastin / self.att_max_me
        self.rec_mod_me = self.lambda_elastin / self.att_mod_me
        self.a_me = self.rec_min_me
        self.b_me = self.rec_max_me
        self.c_me = self.rec_mod_me

        # Adventitia attachment stretch distribution
        self.att_min_ad = 0.8
        self.att_mod_ad = 0.9
        self.att_max_ad = 0.99999

        # Adventitia recruitment stretch distribution
        self.rec_max_ad = self.lambda_elastin / self.att_min_ad
        self.rec_min_ad = self.lambda_elastin / self.att_max_ad
        self.rec_mod_ad = self.lambda_elastin / self.att_mod_ad
        self.a_ad = self.rec_min_ad
        self.b_ad = self.rec_max_ad
        self.c_ad = self.rec_mod_ad

        self.tgf_beta_levels = {
            "TT": 0.713,
            "TC": 0.916,
            "CC": 1.119,
        }  # Genotype-specific latent TGF-Beta production levels (Dr. Mark Bakker at the Utrecht UMC).

        self.smc_mean_fractions = {
            0: 0.7262666,
            1: 0.7132891,
            2: 0.6736118,
            3: 0.7635375,
            4: 0.4898171,
        }  # Mean SMC fractions for different polygenic scores (0-4).

        self.smc_sd_fractions = {
            0: 0.1118832,
            1: 0.1333362,
            2: 0.1479770,
            3: 0.0969533,
            4: 0.1184857,
        }  # Standard deviations of SMC fractions for different polygenic scores (0-4).

        self.polygenic_score = polygenic_score
        self.genotype = genotype

        def _is_missing(v):
            return v is None or (isinstance(v, float) and np.isnan(v))

        # Set TGF-Beta level based on genotype or provided value
        if tgf_beta_level is not None:
            self.tgf_beta_level = tgf_beta_level
        elif _is_missing(genotype):
            self.tgf_beta_level = 1.0
        elif genotype in self.tgf_beta_levels:
            self.tgf_beta_level = self.tgf_beta_levels[genotype]
        else:
            self.tgf_beta_level = 1.0

        # Set SMC fraction based on polygenic score or provided value
        if smc_fraction is not None:
            self.smc_fraction = smc_fraction
        else:
            if _is_missing(polygenic_score):
                score = 0
            else:
                score = int(np.clip(polygenic_score, 0, 4))
            self.smc_fraction = self.smc_mean_fractions[score]

        self.tgf_spike_amount = 0.65  # Amount of TGF-Beta spike due to treatment

        # Immune cell related rates
        self.r_e = 1.0  # Elastin degradation rate by immune cell proteases (years^-1)
        self.r_cm = (
            1.0  # Medial collagen degradation rate by immune cell proteases (years^-1)
        )
        self.t_i0 = 40  # Time for onset of immune cell infiltration (years)
        self.t_treat = 45  # Time for TGF-Beta treatment (years)
        self.i_0 = 0  # Initial level of immune cells in the arterial wall
        self.i_max = 1.0  # Maximum level of immune cells in the arterial wall
        self.k_i = (
            1.25  # Time for immune cell levels to reach half of maximum level (years)
        )
        self.r_pc1 = 1.0  # Collagenase secretion rate by immune cells (years^-1)
        self.r_pc2 = 1.0  # Baseline immune cell collagenase degradation rate (years^-1)
        self.r_pe1 = 1.0  # Elastase secretion rate by immune cells (years^-1)
        self.r_pe2 = 1.0  # Baseline immune cell elastase degradation rate (years^-1)

        # Fibroblast rates
        self.r_f1 = (
            1.0  # Baseline fibroblast migration and proliferation rate (years^-1)
        )
        self.r_f2 = (
            0.5  # Fibroblast population dynamics sensitivity to TGF-Beta (years^-1)
        )
        self.r_f3 = 1.0  # Fibroblast cell death rate (years^-1)

        # Smooth muscle cell rates
        self.beta1_smc = 0.5  # Rate of change in SMCs according to stretch
        self.beta2_smc = 0.0  # Rate of change in SMCs according to change in elastin
        self.beta3_smc = 1.0  # Rate of change in SMCs according to immune cells
        self.beta_wss_smc = 1.0  # or 500
        self.tau_homeo = 1.0
        self.k_active_smc = 12e8  # Material parameter for the active response of SMCs
        self.k_passive_smc = (
            11.8e3  # Material parameter for the passive response of SMCs
        )
        self.lambda_smc_max = 1.4  # Stretch where active force is max
        self.lambda_smc_zero = 2.0  # Stretch limit where active force becomes 0
        self.lambda_att_smc = 1.1  # Attachment stretch of SMCs

        # Procollagen rates
        self.r_p1 = 1.0  # Baseline procollagen secretion rate by fibroblasts (years^-1)
        self.r_p2 = (
            0.5  # Fibroblast procollagen secretion sensitivity to TGF-Beta (years^-1)
        )
        self.r_p3 = 1.0  # Combined baseline procollagen degradation and modification rate (years^-1)

        # Adventitial collagen rates
        self.r_c1 = 0.5  # Baseline adventitial collagen maturation rate (years^-1)
        self.r_c2 = 0.5  # Adventitial collagen degradation rate (years^-1)

        # Zymogen rates
        self.r_z1 = 1.0  # Baseline zymogen secretion rate by fibroblasts (years^-1)
        self.r_z2 = (
            0.5  # Fibroblast zymogen secretion sensitivity to TGF-Beta (years^-1)
        )
        self.r_z3 = 1.0  # Combined baseline zymogen degradation and modification rates (years^-1)

        # Collagenase rates
        self.r_ca1 = 0.5  # Baseline active collagenase maturation rate (years^-1)
        self.r_ca2 = 0.25  # Collagenase inactivation rate by TIMPs (years^-1)
        self.r_ca3 = 0.25  # Inhibitor–collagenase complex formation rate (years^-1)

        # Inhibitor rates
        self.r_i1 = 1.0  # Baseline inhibitor secretion rate (years^-1)
        self.r_i2 = 0.5  # Fibroblast collagenase inhibitor secretion sensitivity to TGF-β (years^-1)
        self.r_i3 = 0.75  # Baseline collagenase inhibitor degradation rate (years^-1)
        self.r_i4 = 0.25  # Inhibitor–collagenase complex formation rate (duplicate of r_ca3) (years^-1)

        # TGF-Beta rates
        self.r_betal1 = (
            0.1  # Fibroblast latent TGF-Beta secretion sensitivity to active TGF-Beta
        )
        self.r_betal2 = 5.0  # Fibroblast latent TGF-Beta secretion sensitivity to deviations from mechanical homeostasis (Parameter study with [0.1, 1.0, 5.0, 10.0])
        self.r_betal3 = (
            1.0  # Fibroblast latent TGF-Beta secretion sensitivity to collagen levels
        )
        self.r_betal4 = 1.0  # Combined baseline latent TGF-Bet degradation/modification rate (years^-1)
        self.r_betal5 = 1.0  # Latent TGF-Beta modification rate by integrin/ECM/Stretch–dependent mechanism (years^-1)
        self.r_beta1 = 0.5  # Baseline latent TGF-Beta activation rate
        self.r_beta2 = 1.0  # Modification rate by integrin/ECM/Stretch–dependent mechanism (same as r_betaL5)
        self.r_beta3 = 1.0  # Baseline active TGF-Beta degradation rate (years^-1)

        # Mechanical model parameters
        self.remodel_time = (
            10  # Averaging time period for attachment stretch remodelling
        )
        self.t_sim = 90  # Simulation time in years
        self.width_att_dist = (
            self.att_max_ad - self.att_min_ad
        )  # Width of the attachment stretch distribution
        self.skew_att_dist = (
            self.att_mod_ad - self.att_min_ad
        ) / self.width_att_dist  # Skew of the attachment stretch distribution
        self.width_att_dist_me = (
            self.att_max_me - self.att_min_me
        )  # Width of the attachment stretch distribution for the media
        self.skew_att_dist_me = (
            (self.att_mod_me - self.att_min_me) / self.width_att_dist_me
        )  # Skew of the attachment stretch distribution for the media

        # Set initial values for densities of cells and proteins in the arterial wall
        self.alpha_init = 0.1
        self.init_fibroblast = 1.0
        self.init_muscle_cells = 1.0
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
        self.init_latent_tgf_beta = 0.0
        self.init_active_tgf_beta = 0.0

        self.refresh_physics()

    def refresh_physics(self):
        """
        Recomputes load-bearing fractions and material stiffnesses (collagen, elastin, vSMCs).

        This method refreshes the media and adventitia recruitment-stretch distributions from
        the fixed attachment stretches, partitions the load-bearing fractions of the arterial wall constituents,
        and calibrates the stiffness parameters (k_collagen, k_elastin, k_muscle_p, k_muscle_a) so that each constituent's
        Cauchy stress carries its assigned share of the systolic pressure via the thin-walled (Laplace) force balance.
        """
        self.att_max_me_current = self.att_max_me

        # Refresh media recruitment-stretch distribution vertices
        self.rec_max_me = self.lambda_elastin / self.att_min_me
        self.rec_min_me = self.lambda_elastin / self.att_max_me
        self.rec_mod_me = self.lambda_elastin / self.att_mod_me
        self.a_me, self.b_me, self.c_me = (
            self.rec_min_me,
            self.rec_max_me,
            self.rec_mod_me,
        )

        # Refresh adventitia recruitment-stretch distribution vertices
        self.rec_max_ad = self.lambda_elastin / self.att_min_ad
        self.rec_min_ad = self.lambda_elastin / self.att_max_ad
        self.rec_mod_ad = self.lambda_elastin / self.att_mod_ad
        self.a_ad, self.b_ad, self.c_ad = (
            self.rec_min_ad,
            self.rec_max_ad,
            self.rec_mod_ad,
        )

        # Partition load-bearing fractions of the arterial wall constituents
        # SMC contributes smc_fraction/2 each to passive and active response;
        # the remaining splits 1:2 between elastin and collagen based in literature.
        self.load_borne_muscle_p = self.smc_fraction / 2
        self.load_borne_muscle_a = self.load_borne_muscle_p
        self.load_borne_elastin = (1 / 3) * (
            1 - (self.load_borne_muscle_p + self.load_borne_muscle_a)
        )
        self.load_borne_collagen = 1.0 - (
            self.load_borne_elastin
            + self.load_borne_muscle_p
            + self.load_borne_muscle_a
        )

        # Common Laplace force-balance factor (P * R0 * lambda^2 * lambda_z / H)
        # per layer, used to convert load-bearing fractions into stiffnesses.
        self.common_factor_me = (
            self.pressure_sys
            * self.radius_tzero
            * self.lambda_elastin**2
            * self.lambda_z
            / self.thickness_me
        )
        self.common_factor_ad = (
            self.pressure_sys
            * self.radius_tzero
            * self.lambda_elastin**2
            * self.lambda_z
            / self.thickness_ad
        )

        # Integrate the piecewise triangular collagen recruitment distribution
        # at the homeostatic stretch x to obtain the collagen stress denominator.
        x = self.lambda_elastin
        v_a = self.a_me
        v_b = self.b_me
        v_c = self.c_me

        gamma_unit = 1.0 / ((v_b - v_a) * (v_c - v_a))
        delta_unit = 1.0 / ((v_b - v_a) * (v_b - v_c))

        if x < v_a:
            # Below recruitment: no collagen load bearing (guard against /0).
            collagen_denominator = 1e-9
        elif x < v_c:
            # Min -> mode branch of the triangular distribution.
            collagen_denominator = (
                x * gamma_unit * 2 * ((x + v_a) * np.log(x / v_a) + 2 * (v_a - x))
            )
        elif x <= v_b:
            # Mode -> max branch of the triangular distribution.
            term1 = (x + v_a) * np.log(v_c / v_a) + v_a - v_c + ((v_a - v_c) / v_c) * x
            term2 = (x + v_b) * np.log(x / v_c) + v_b + v_c - ((v_b + v_c) / v_c) * x
            collagen_denominator = (
                x * gamma_unit * 2 * term1 - x * delta_unit * 2 * term2
            )
        else:
            # Above max recruitment, all collagen is load bearing.
            term1 = (x + v_a) * np.log(v_c / v_a) + v_a - v_c + ((v_a - v_c) / v_c) * x
            term2 = (x + v_b) * np.log(v_b / v_c) - v_b + v_c - ((v_b - v_c) / v_c) * x
            collagen_denominator = (
                x * gamma_unit * 2 * term1 - x * delta_unit * 2 * term2
            )

        # Collagen stiffness parameter
        self.k_collagen = (
            self.load_borne_collagen * (self.common_factor_me) / collagen_denominator
        )

        # Triangular-distribution stress coefficients per layer
        self.gamma_me = self.k_collagen / (
            (self.b_me - self.a_me) * (self.c_me - self.a_me)
        )
        self.delta_me = self.k_collagen / (
            (self.b_me - self.a_me) * (self.b_me - self.c_me)
        )
        self.gamma_ad = (
            self.k_collagen
            * self.collagen_ratio_ad_me
            / ((self.b_ad - self.a_ad) * (self.c_ad - self.a_ad))
        )
        self.delta_ad = (
            self.k_collagen
            * self.collagen_ratio_ad_me
            / ((self.b_ad - self.a_ad) * (self.b_ad - self.c_ad))
        )

        # Elastin (neo-Hookean) stiffness parameter
        self.k_elastin = (
            self.load_borne_elastin
            * (self.common_factor_me)
            / (
                self.lambda_elastin**2
                * (1 - (1 / (self.lambda_z**2 * self.lambda_elastin**4)))
            )
        )

        muscle_a_denominator = (
            self.vasodil_conc
            * self.lambda_muscle
            * (
                1
                - (
                    (self.musc_mean - self.lambda_muscle)
                    / (self.musc_mean - self.musc_min)
                )
                ** 2
            )
        )

        # Active SMC (length-tension) stiffness parameter
        self.k_muscle_a = (
            self.load_borne_muscle_a * self.common_factor_me / muscle_a_denominator
        )

        # Passive SMC (neo-Hookean) stiffness parameter
        self.k_muscle_p = (
            self.load_borne_muscle_p
            * (self.common_factor_me)
            / (
                self.lambda_muscle**2
                * (1 - 1 / (self.lambda_z**2 * self.lambda_muscle**4))
            )
        )

        print(
            f"Updated physics: K elastin: {self.k_elastin:.2f}, K collagen: {self.k_collagen:.2f}, K muscle passive: {self.k_muscle_p:.2f}, K muscle active: {self.k_muscle_a:.2f}"
        )

    def to_dict(self):
        return self.__dict__
