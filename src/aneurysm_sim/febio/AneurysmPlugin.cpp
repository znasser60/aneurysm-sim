#include <FEBioMech/FEUncoupledMaterial.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FECore/FECoreKernel.h>
#include <FECore/log.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstring>


static double MaxPrincipalStretch(const mat3d& F)
{
    /**
     * @brief  Calculates maximum principal stretch:
     * Uses the deformation tensor to calculate the 
     * maximum principal stretch through eigen-decomposition
     * of C = F^TF. 
     * 
     * @param  F:  Deformation tensor
     * 
     * @return The maximum principal stretch
     */
    mat3ds C = (F.transpose() * F).sym();
    double eval[3];
    C.eigen(eval);
    double lam_sq_max = std::max({eval[0], eval[1], eval[2]});
    return sqrt(lam_sq_max);
}


class FEAneurysmPoint :
public FEMaterialPointData
{
public:
    double bio_time;
    double bio_last_update_time;
    double homeostatic_stretch;

    double immune_cells;
    double elastases_imm;
    double collagenases_imm;
    double elastin_me;
    double collagen_me;
    double muscle_cells;
    double fibroblast;
    double procollagen;
    double collagenase_fib;
    double zymogen;
    double timp;
    double latent_tgf_beta;
    double active_tgf_beta;
    double collagen_ad;

    double lambda_rec_min;
    double lambda_rec_max;
    double lambda_rec_mode;

    static const int MAX_HIST = 2000;
    double lc_max_history[MAX_HIST];
    int    hist_count;

    FEAneurysmPoint(FEMaterialPointData* pt) : FEMaterialPointData(pt) {}

    FEMaterialPointData* Copy() override {
        FEAneurysmPoint* p = new FEAneurysmPoint(*this);
        if (m_pNext) p->m_pNext = m_pNext->Copy();
        return p;
    }

    void Serialize(DumpStream& ar) override {
        FEMaterialPointData::Serialize(ar);
        ar & bio_time & bio_last_update_time & homeostatic_stretch;
        ar & immune_cells & elastases_imm & collagenases_imm;
        ar & elastin_me & collagen_me & muscle_cells;
        ar & fibroblast & procollagen & collagenase_fib & zymogen & timp;
        ar & latent_tgf_beta & active_tgf_beta & collagen_ad;
        ar & lambda_rec_min & lambda_rec_max & lambda_rec_mode;
        ar & hist_count;
        for (int i = 0; i < MAX_HIST; ++i) ar & lc_max_history[i];
    }

    void Init() override {
        bio_time             = 0.0;
        bio_last_update_time = -1.0;
        homeostatic_stretch  = -1.0;   // sentinel: negative = not yet initialised

        immune_cells      = 0.0;
        elastases_imm     = 0.0;
        collagenases_imm  = 0.0;
        elastin_me        = 1.0;
        collagen_me       = 1.0;
        muscle_cells      = 1.0;

        fibroblast        = 1.0;
        procollagen       = 1.0;
        zymogen           = 1.0;
        timp              = 1.0;
        collagenase_fib   = 1.0;
        latent_tgf_beta   = 0.0;
        active_tgf_beta   = 0.0;
        collagen_ad       = 1.0;

        lambda_rec_min  = 1.3 / 1.0;  
        lambda_rec_max  = 1.3 / 0.9; 
        lambda_rec_mode = 1.3 / 0.95; 

        hist_count = 0;
        for (int i = 0; i < MAX_HIST; ++i) lc_max_history[i] = 1.3;

        FEMaterialPointData::Init();
    }
};

const int FEAneurysmPoint::MAX_HIST;


class FEAneurysmMaterial :
public FEUncoupledMaterial
{
public:
    FEUncoupledMaterial* m_pPassive;
    FEUncoupledMaterial* m_pActive;

    // Layer and genotype identifiers
    int    layer; // 0 = Media, 1 = Adventitia
    int    genotype; // 0 = TT, 1 = TC, 2 = CC
    int    polygenic_score; // 0–4
    double tgf_beta_level;
    double smc_fraction;

    // Simulation timing
    double febio_disease_start;
    double years_per_febio_sec;

    // Immune infiltration patch geometry
    double patch_x, patch_y, patch_z;
    double spread_sigma;

    // Immune / protease
    double r_e    = 1.0;
    double r_cm   = 1.0;
    double t_i0   = 40.0;
    double i_0    = 0.0;
    double i_max  = 1.0;
    double k_i    = 1.25;
    double r_pc1  = 1.0;
    double r_pc2  = 1.0;

    // Fibroblast
    double r_f1 = 1.0;
    double r_f2 = 0.5;
    double r_f3 = 1.0;

    // Procollagen
    double r_p1 = 1.0;
    double r_p2 = 0.5;
    double r_p3 = 1.0;

    // Adventitial collagen maturation/degradation
    double r_c1 = 0.5;
    double r_c2 = 0.5;

    // Zymogen
    double r_z1 = 1.0;
    double r_z2 = 0.5;
    double r_z3 = 1.0;

    // Collagenase
    double r_ca1 = 0.5;
    double r_ca2 = 0.25;
    double r_ca3 = 0.25;

    // TIMP / inhibitor
    double r_i1 = 1.0;
    double r_i2 = 0.5;
    double r_i3 = 0.75;
    double r_i4 = 0.25;

    // TGF-beta latent
    double r_betal1 = 0.1;
    double r_betal2 = 5.0;
    double r_betal3 = 1.0;
    double r_betal4 = 1.0;
    double r_betal5 = 1.0;

    // TGF-beta active
    double r_beta1 = 0.5;
    double r_beta2 = 1.0;
    double r_beta3 = 1.0;

    // SMC
    double beta1_smc     = 0.5;
    double beta2_smc     = 0.0;
    double beta3_smc     = 1.0;
    double lambda_att_smc = 1.1;
    // double lam_sys = 1.3; 
    // double lam_muscle_init = 1.13;
    // double rec_muscle = lam_sys / lam_muscle_init;

    // Remodelling parameters
    double remodel_time    = 10.0;  // years (for reference only, N is fixed below)
    int    remodel_N       = 1000;  // fixed step count matching 1D dt=0.01 yr window. Equal to 10 year remodel.

    // Adventitia attachment stretch distribution
    double width_att_dist  = 0.19999;  // att_max_ad - att_min_ad = 0.99999 - 0.8
    double skew_att_dist   = 0.500;    // (att_mod_ad - att_min_ad) / width

    // Media attachment stretch distribution 
    double width_att_dist_me = 0.06999; // att_max_me - att_min_me = 1.07 - 1.00001
    double skew_att_dist_me  = 0.143;   // (att_mod_me - att_min_me) / width_me

    // Remodelling rate
    // media uses constant alpha_CM0 
    double alpha_init = 0.1; 

    // Sub-stepping resolution to match 1D model resolution
    double dt_ref = 0.01;   // yr

    double febio_treatment_start = -1.0;  // negative = disabled
    double tgf_spike_amount      = 0.65;        


    FEAneurysmMaterial(FEModel* fem) : FEUncoupledMaterial(fem)
    {
        m_pPassive = nullptr;
        m_pActive  = nullptr;

        layer           = 0;
        genotype        = 0;
        polygenic_score = 0;

        febio_disease_start = 2.0;
        years_per_febio_sec = 45.0;

        patch_x = 1.27; patch_y = 0.0; patch_z = 5.0;
        spread_sigma = 0.3;
    }

    DECLARE_FECORE_CLASS();

    bool Init() override
    {
        const double tgf_map[3] = {0.713, 0.916, 1.119};
        const double smc_map[5] = {0.7262666, 0.7132891, 0.6736118, 0.7635375, 0.4898171};

        genotype        = std::max(0, std::min(genotype, 2));
        polygenic_score = std::max(0, std::min(polygenic_score, 4));

        tgf_beta_level = tgf_map[genotype];
        smc_fraction   = smc_map[polygenic_score];

        return FEUncoupledMaterial::Init();
    }

    FEMaterialPointData* CreateMaterialPointData() override
    {
        FEMaterialPointData* pd = m_pPassive->CreateMaterialPointData();
        if (m_pActive) pd->Append(m_pActive->CreateMaterialPointData());
        return new FEAneurysmPoint(pd);
    }

    double ComputeScale(const FEAneurysmPoint& pt) const
    {
        if (layer == 0) {
            double w_e  = (1.0 / 3.0) * (1.0 - smc_fraction);
            double w_cm = (2.0 / 3.0) * (1.0 - smc_fraction);
            double w_mp = 0.5 * smc_fraction;
            return (w_e * pt.elastin_me + w_cm * pt.collagen_me + w_mp * pt.muscle_cells)
                 / (w_e + w_cm + w_mp);
        }
        else {
            return pt.collagen_ad;
        }
    }

    double GetCircumferentialStretch(FEMaterialPoint& mp) const
    {
        FEElasticMaterialPoint* ep = mp.ExtractData<FEElasticMaterialPoint>();
        if (!ep) return 1.0;
        return MaxPrincipalStretch(ep->m_F);
    }

    double CalcLamAttMax(FEAneurysmPoint& pt) const
    {
        int N     = std::min(remodel_N, FEAneurysmPoint::MAX_HIST);
        int count = pt.hist_count;

        if (count == 0)
            return (pt.lambda_rec_min > 0) ? 1.3 / pt.lambda_rec_min : 1.0;

        double baseline = pt.lc_max_history[0];  // homeostatic lc_max, set pre-disease

        if (count < N) {
            int    missing   = N - count;
            double sum_ghost = missing * baseline;
            double sum_rec   = 0.0;
            for (int k = 0; k < count; ++k) sum_rec += pt.lc_max_history[k];
            return (sum_ghost + sum_rec) / (double)N;
        }
        else {
            int    use   = std::min(count, FEAneurysmPoint::MAX_HIST);
            int    start = use - N;
            double sum   = 0.0;
            for (int k = start; k < use; ++k) sum += pt.lc_max_history[k];
            return sum / (double)N;
        }
    }

    double GetScale(FEMaterialPoint& mp) const
    {
        FEAneurysmPoint& pt = *mp.ExtractData<FEAneurysmPoint>();
        double t_febio = GetTimeInfo().currentTime;
        if (t_febio <= febio_disease_start) return 1.0;
        if (pt.bio_last_update_time < 0)    return 1.0;
        return ComputeScale(pt);
    }

    mat3ds DevStress(FEMaterialPoint& mp) override
    {
        if (!m_pPassive) return mat3ds(0.0);

        FEAneurysmPoint& pt = *mp.ExtractData<FEAneurysmPoint>();
        double t_febio  = GetTimeInfo().currentTime;
        double dt_febio = GetTimeInfo().timeIncrement;

        double lam = GetCircumferentialStretch(mp);

        mat3ds s_pass_raw = m_pPassive->DevStress(mp);
        mat3ds s_act_raw(0.0);
        if (m_pActive) s_act_raw = m_pActive->DevStress(mp);

        if (t_febio <= febio_disease_start) {

            if (pt.homeostatic_stretch < 0.0) {
                if (layer == 0) {
                    pt.lambda_rec_min  = 1.3 / 1.07;       // 1.21495
                    pt.lambda_rec_max  = 1.3 / 1.00001;    // 1.29999
                    pt.lambda_rec_mode = 1.3 / 1.01;        // 1.28713
                }
            }

            pt.homeostatic_stretch  = lam;
            pt.lc_max_history[0]    = lam / pt.lambda_rec_min;
            pt.hist_count           = 1;
            return s_pass_raw + s_act_raw;
        }

        if (pt.bio_last_update_time < 0) {
            pt.bio_last_update_time = t_febio;
            return s_pass_raw + s_act_raw;
        }
        if (std::fabs(t_febio - pt.bio_last_update_time) > 0) {
            pt.bio_last_update_time = t_febio;

            double dt_bio = dt_febio * years_per_febio_sec;
            pt.bio_time  += dt_bio;

            vec3d  r0      = mp.m_r0;
            double dist_sq = pow(r0.x - patch_x, 2)
                           + pow(r0.y - patch_y, 2)
                           + pow(r0.z - patch_z, 2);
            double gw         = exp(-dist_sq / (2.0 * spread_sigma * spread_sigma));
            double local_imax = i_max * gw;

            if (pt.bio_time >= t_i0) {
                pt.immune_cells = i_0 + local_imax
                    * ((pt.bio_time - t_i0) / (k_i + (pt.bio_time - t_i0)));
            }

            double lc_max  = lam / pt.lambda_rec_min;
            double lc_min  = lam / pt.lambda_rec_max;
            double lc_mode = lam / pt.lambda_rec_mode;

            double lam_att_max = CalcLamAttMax(pt);

            double width = (layer == 1) ? width_att_dist    : width_att_dist_me;
            double skew  = (layer == 1) ? skew_att_dist     : skew_att_dist_me;

            double lam_att_min  = lam_att_max - width;
            double lam_att_mode = lam_att_min + skew * (lam_att_max - lam_att_min);

            double f_lam = (lc_max <= lam_att_max)
                         ? 0.0
                         : (lc_max - lam_att_max) / lam_att_max;

            double alpha_adv = alpha_init * (pt.fibroblast / pt.collagen_ad)
                             * sqrt(pt.collagen_ad * pt.collagenase_fib);
            double alpha_med = alpha_init;
            // double lam_muscle = lam / rec_muscle;
            double eps_str  = std::max(0.0, (lam - pt.homeostatic_stretch)
                              / pt.homeostatic_stretch);
            double eps_elas = pt.elastin_me - 1.0;
            double eps_imm  = i_0 - pt.immune_cells;

            auto safe_div = [](double num, double den) -> double {
                return (std::fabs(den) > 1e-12) ? num / den : 0.0;
            };

            // if (dist_sq < 0.002 && std::fmod(pt.bio_time, 1.0) < dt_bio) {
            //     feLog("Geno:%d layer:%d t:%.1f lam:%.4f lc_max:%.4f att:%.4f "
            //           "f_lam:%.5f col_ad:%.4f col_me:%.4f elastin:%.4f "
            //           "muscle:%.4f lat_tgf:%.5f act_tgf:%.5f "
            //           "sxx:%.4f syy:%.4f szz:%.4f scale:%.4f\n",
            //         genotype, layer, pt.bio_time, lam, lc_max, lam_att_max,
            //         f_lam, pt.collagen_ad, pt.collagen_me, pt.elastin_me,
            //         pt.muscle_cells, pt.latent_tgf_beta, pt.active_tgf_beta,
            //         s_pass_raw.xx(), s_pass_raw.yy(), s_pass_raw.zz(),
            //         ComputeScale(pt));
            // }

            int    n_sub  = std::max(1, (int)std::round(dt_bio / dt_ref));
            double dt_sub = dt_bio / n_sub;

            for (int s = 0; s < n_sub; ++s) {
                // snapshot
                double elas_0    = pt.elastases_imm;
                double col_imm_0 = pt.collagenases_imm;
                double elast_0   = pt.elastin_me;
                double colme_0   = pt.collagen_me;
                double smc_0     = pt.muscle_cells;
                double fib_0     = pt.fibroblast;
                double procol_0  = pt.procollagen;
                double zym_0     = pt.zymogen;
                double ca_0      = pt.collagenase_fib;
                double timp_0    = pt.timp;
                double colad_0   = pt.collagen_ad;
                double lat_0     = pt.latent_tgf_beta;
                double act_0     = pt.active_tgf_beta;

                pt.elastases_imm    += dt_sub * (r_pc1 * pt.immune_cells - r_pc2 * elas_0);
                pt.collagenases_imm += dt_sub * (r_pc1 * pt.immune_cells - r_pc2 * col_imm_0);
                pt.elastin_me       *= exp(-r_e  * elas_0  * dt_sub);
                pt.collagen_me      *= exp(-r_cm * col_imm_0 * dt_sub);
                pt.muscle_cells     += dt_sub * smc_0
                                    * (beta1_smc * eps_str + beta2_smc * eps_elas + beta3_smc * eps_imm);
                pt.fibroblast       += dt_sub * ((r_f1 + r_f2 * act_0) * fib_0 - r_f3 * fib_0);
                pt.procollagen      += dt_sub * ((r_p1 + r_p2 * act_0) * fib_0 - r_p3 * procol_0);
                pt.zymogen          += dt_sub * ((r_z1 / (1.0 + r_z2 * act_0)) * fib_0 - r_z3 * zym_0);
                pt.collagenase_fib  += dt_sub * (r_ca1 * zym_0 - (r_ca2 + r_ca3 * timp_0) * ca_0);
                pt.timp             += dt_sub * ((r_i1 + r_i2 * act_0) * fib_0
                                                - (r_i3 + r_i4 * ca_0) * timp_0);
                pt.collagen_ad      += dt_sub * (r_c1 * procol_0 - r_c2 * ca_0 * colad_0);

                double term1 = (r_betal1 * act_0 + r_betal2 * f_lam) / (1.0 + r_betal3 * colad_0);
                double term3 = r_betal4 + r_betal5 * f_lam * fib_0;
                pt.latent_tgf_beta  += dt_sub * (term1 * fib_0 * tgf_beta_level - term3 * lat_0);
                pt.active_tgf_beta  += dt_sub * ((r_beta1 + r_beta2 * f_lam * fib_0) * lat_0
                                                - r_beta3 * act_0);

                
                // TGF-beta treatment spike
                if (febio_treatment_start > 0.0) {
                    double t_treat_bio = (febio_treatment_start - febio_disease_start) 
                                        * years_per_febio_sec;
                    if (pt.bio_time >= t_treat_bio && 
                        pt.bio_time - dt_sub < t_treat_bio) {
                        pt.active_tgf_beta += tgf_spike_amount;
                    }
                }

                if (layer == 1) {
                    pt.lambda_rec_min  += dt_sub * alpha_adv
                        * safe_div(lc_max - lam_att_max, lam_att_max);
                    pt.lambda_rec_max  += dt_sub * alpha_adv
                        * safe_div(lc_min  - lam_att_min,  lam_att_min);
                    pt.lambda_rec_mode += dt_sub * alpha_adv
                        * safe_div(lc_mode - lam_att_mode, lam_att_mode);
                }
                else {
                    pt.lambda_rec_min  += dt_sub * alpha_med
                        * safe_div(lc_max - lam_att_max, lam_att_max);
                    pt.lambda_rec_max  += dt_sub * alpha_med
                        * safe_div(lc_min  - lam_att_min,  lam_att_min);
                    pt.lambda_rec_mode += dt_sub * alpha_med
                        * safe_div(lc_mode - lam_att_mode, lam_att_mode);
                }

            } // end sub-loop

            // --- Append lc_max to history (once per FEBio step, after sub-loop) ---
            if (pt.hist_count < FEAneurysmPoint::MAX_HIST) {
                pt.lc_max_history[pt.hist_count] = lc_max;
            }
            else {
                std::memmove(&pt.lc_max_history[0], &pt.lc_max_history[1],
                             (FEAneurysmPoint::MAX_HIST - 1) * sizeof(double));
                pt.lc_max_history[FEAneurysmPoint::MAX_HIST - 1] = lc_max;
            }
            pt.hist_count++;

        } 

        double passive_scale = ComputeScale(pt);
        if (layer == 0) {
            return (s_pass_raw * passive_scale) + (s_act_raw * pt.muscle_cells);
        }
        else {
            return s_pass_raw * passive_scale;
        }
    }

    tens4ds DevTangent(FEMaterialPoint& mp) override
    {
        if (!m_pPassive) return tens4ds(0.0);

        FEAneurysmPoint& pt = *mp.ExtractData<FEAneurysmPoint>();
        double t_febio = GetTimeInfo().currentTime;

        tens4ds c_pass_raw = m_pPassive->DevTangent(mp);
        tens4ds c_act_raw(0.0);
        if (m_pActive) c_act_raw = m_pActive->DevTangent(mp);

        if (t_febio <= febio_disease_start)
            return c_pass_raw + c_act_raw;

        double passive_scale = ComputeScale(pt);
        if (layer == 0) {
            return (c_pass_raw * passive_scale) + (c_act_raw * pt.muscle_cells);
        }
        else {
            return c_pass_raw * passive_scale;
        }
    }
};

BEGIN_FECORE_CLASS(FEAneurysmMaterial, FEUncoupledMaterial)
    ADD_PROPERTY(m_pPassive, "passive_solid");
    ADD_PROPERTY(m_pActive,  "active_solid", FEProperty::Optional);
    ADD_PARAMETER(layer,               "layer");
    ADD_PARAMETER(genotype,            "genotype");
    ADD_PARAMETER(polygenic_score,     "polygenic_score");
    ADD_PARAMETER(febio_disease_start, "febio_disease_start");
    ADD_PARAMETER(years_per_febio_sec, "years_per_febio_sec");
    ADD_PARAMETER(patch_x,             "patch_x");
    ADD_PARAMETER(patch_y,             "patch_y");
    ADD_PARAMETER(patch_z,             "patch_z");
    ADD_PARAMETER(spread_sigma,        "spread_sigma");
    ADD_PARAMETER(remodel_N,           "remodel_N");
    ADD_PARAMETER(r_betal2,            "r_betal2");
    ADD_PARAMETER(febio_treatment_start, "febio_treatment_start");
    ADD_PARAMETER(tgf_spike_amount,      "tgf_spike_amount");
END_FECORE_CLASS();


extern "C" {

    FECORE_EXPORT int GetSDKVersion() {
        return FE_SDK_VERSION;
    }

    FECORE_EXPORT bool PluginInitialize(FECoreKernel& fecore) {
        REGISTER_FECORE_CLASS(FEAneurysmMaterial, "aneurysm_0D_coupled");
        return true;
    }

}