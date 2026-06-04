#include <FEBioMech/FEUncoupledMaterial.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FECore/FECoreKernel.h>
#include <cmath>
#include <algorithm>
#include <vector>


static double MaxPrincipalStretch(const mat3d& F)
{
    mat3ds C = (F.transpose() * F).sym();
    double eval[3];
    C.eigen(eval);  
    double lam_sq_max = std::max({eval[0], eval[1], eval[2]});
    double lam_max = sqrt(lam_sq_max);
    return lam_max;
}


class FEAneurysmPoint: 
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

    static const int MAX_HIST = 500;
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
        bio_time = 0.0;
        bio_last_update_time = -1.0;
        homeostatic_stretch = -1.0;

        immune_cells = 0.0;
        elastases_imm = 0.0;
        collagenases_imm = 0.0;
        elastin_me = 1.0;
        collagen_me = 1.0;
        muscle_cells = 1.0;

        fibroblast = 1.0;
        procollagen = 1.0;
        zymogen = 1.0;
        timp = 1.0;
        latent_tgf_beta = 0.0;
        active_tgf_beta = 0.0;
        collagen_ad = 1.0;

        collagenase_fib = 1.0;  

        // lambda_rec = c_lambda_elastin / c_att
        lambda_rec_min = 1.3 / 0.99999;
        lambda_rec_max = 1.3 / 0.8;
        lambda_rec_mode = 1.3 / 0.9;

        hist_count = 0;
        for (int i = 0; i < MAX_HIST; ++i) lc_max_history[i] = 1.3;

        FEMaterialPointData::Init();
    }
};

const int FEAneurysmPoint::MAX_HIST;
class FEAneurysmMaterial : public FEUncoupledMaterial
{
public:
    // FEUncoupledMaterial* m_base_mat;
    FEUncoupledMaterial* m_pPassive;
    FEUncoupledMaterial* m_pActive;

    int    layer; // 0 = Media, 1 = Adventitia
    int    genotype; // 0 = TT, 1 = TC, 2 = CC
    int    polygenic_score; // 0–4

    double tgf_beta_level;
    double smc_fraction;

    double febio_disease_start;
    double years_per_febio_sec;

    double patch_x, patch_y, patch_z;
    double spread_sigma;

    double r_e    = 1.0,  r_cm   = 1.0;
    double t_i0   = 40.0, i_0    = 0.0, i_max = 1.0, k_i = 1.25;
    double r_pc1  = 1.0,  r_pc2  = 1.0;
    double r_f1   = 1.0,  r_f2   = 0.5,  r_f3   = 1.0;
    double r_p1   = 1.0,  r_p2   = 0.5,  r_p3   = 1.0;
    double r_c1   = 0.5,  r_c2   = 0.5;
    double r_z1   = 1.0,  r_z2   = 0.5,  r_z3   = 1.0;
    double r_ca1  = 0.5,  r_ca2  = 0.25, r_ca3  = 0.25;
    double r_i1   = 1.0,  r_i2   = 0.5,  r_i3   = 0.75, r_i4 = 0.25;
    double r_betal1 = 0.1, r_betal2 = 5.0, r_betal3 = 1.0;
    double r_betal4 = 1.0, r_betal5 = 1.0;
    double r_beta1  = 0.5, r_beta2  = 1.0, r_beta3 = 1.0;
    double beta1_smc = 0.75, beta2_smc = 0.0, beta3_smc = 1.0;
    double lambda_att_smc   = 1.1;

    double remodel_time   = 10.0;
    double width_att_dist = 0.1;
    double skew_att_dist  = 0.5;
    double alpha_init     = 1.15;

    FEAneurysmMaterial(FEModel* fem) : FEUncoupledMaterial(fem)
    {
        m_pPassive          = nullptr;
        m_pActive           = nullptr;
        layer               = 0;
        genotype            = 0;
        polygenic_score     = 0;
        febio_disease_start = 2.0;
        years_per_febio_sec = 45.0;
        patch_x = 1.27; patch_y = 0.0; patch_z = 5.0;
        spread_sigma = 0.3;
    }

    DECLARE_FECORE_CLASS();

    bool Init() override
    {
        const double tgf_map[3]  = {0.713, 0.916, 1.119};
        const double smc_map[5]  = {0.7262666, 0.7132891, 0.6736118, 0.7635375, 0.4898171};

        genotype        = std::max(0, std::min(genotype, 2));
        polygenic_score = std::max(0, std::min(polygenic_score, 4));

        tgf_beta_level = tgf_map[genotype];
        smc_fraction   = smc_map[polygenic_score];

        return FEUncoupledMaterial::Init();
    }

    FEMaterialPointData* CreateMaterialPointData() override
    {
        FEMaterialPointData* pd = m_pPassive->CreateMaterialPointData();
        if (m_pActive) {
            pd->Append(m_pActive->CreateMaterialPointData());
        }
        return new FEAneurysmPoint(pd);
    }

    double ComputeScale(const FEAneurysmPoint& pt) const
    {
        if (layer == 0) { 
            double w_e  = (1.0 / 3.0) * (1.0 - smc_fraction);
            double w_cm = (2.0 / 3.0) * (1.0 - smc_fraction);
            double w_mp  = 0.5*smc_fraction;
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

    double CalcLamAttMax(FEAneurysmPoint& pt, double dt_bio) const
    {
        int N = std::max(1, (int)(remodel_time / dt_bio));
        int count = pt.hist_count;
        int use = std::min(count, FEAneurysmPoint::MAX_HIST);
        int start = std::max(0, use - N);
        double sum = 0.0;
        int n_used = 0;
        for (int k = start; k < use; ++k) { sum += pt.lc_max_history[k]; ++n_used; }
        if (n_used == 0) return pt.lambda_rec_min > 0 ? 1.3 / pt.lambda_rec_min : 1.0;
        return (sum * dt_bio) / remodel_time;
    }

    double GetScale(FEMaterialPoint& mp) const
    {
        FEAneurysmPoint& pt = *mp.ExtractData<FEAneurysmPoint>();
        double t_febio = GetTimeInfo().currentTime;
        if (t_febio <= febio_disease_start) return 1.0;
        if (pt.bio_last_update_time < 0) return 1.0;
        return ComputeScale(pt);
    }

    mat3ds DevStress(FEMaterialPoint& mp) override
    {
        if (!m_pPassive) return mat3ds(0.0);

        FEAneurysmPoint& pt = *mp.ExtractData<FEAneurysmPoint>();
        double t_febio = GetTimeInfo().currentTime;
        double dt_febio = GetTimeInfo().timeIncrement;
        double lam = GetCircumferentialStretch(mp); 
        mat3ds s_pass_raw = m_pPassive->DevStress(mp);
        mat3ds s_act_raw(0.0);
        if (m_pActive) s_act_raw = m_pActive->DevStress(mp);

        if (t_febio <= febio_disease_start) {
            pt.homeostatic_stretch = lam;
            for (int i = 0; i < FEAneurysmPoint::MAX_HIST; ++i) pt.lc_max_history[i] = lam / pt.lambda_rec_min;
            return s_pass_raw + s_act_raw;
        }

        if (pt.bio_last_update_time < 0) {
            pt.bio_last_update_time = t_febio;
            return s_pass_raw + s_act_raw;
        }

        if (std::fabs(t_febio - pt.bio_last_update_time) > 0){
            pt.bio_last_update_time = t_febio;

            double dt_bio = dt_febio * years_per_febio_sec;
            pt.bio_time += dt_bio;

            vec3d  r0       = mp.m_r0;
            double dist_sq  = pow(r0.x - patch_x, 2)
                            + pow(r0.y - patch_y, 2)
                            + pow(r0.z - patch_z, 2);
            double gw       = exp(-dist_sq / (2.0 * spread_sigma * spread_sigma));
            double local_imax = i_max * gw;

            if (pt.bio_time >= t_i0) {
                pt.immune_cells = i_0 + local_imax
                    * ((pt.bio_time - t_i0) / (k_i + (pt.bio_time - t_i0)));
            }
            double lam_ref = (pt.homeostatic_stretch > 0.5)
                           ? pt.homeostatic_stretch : lam;
            double f_lam   = std::max(0.0, (lam - lam_ref) / lam_ref);

            pt.elastases_imm += dt_bio * (r_pc1 * pt.immune_cells - r_pc2 * pt.elastases_imm);
            pt.collagenases_imm += dt_bio * (r_pc1 * pt.immune_cells - r_pc2 * pt.collagenases_imm);
            pt.elastin_me *= exp(-r_e * pt.elastases_imm * dt_bio);
            pt.collagen_me *= exp(-r_cm * pt.collagenases_imm * dt_bio);

            // double eps_str  = std::max(0.0, (lam - lambda_att_smc) / lambda_att_smc);
            double eps_str = std::max(0.0, (lam - pt.homeostatic_stretch) / pt.homeostatic_stretch);
            double eps_elas = (pt.elastin_me  - 1.0);
            double eps_imm  = (i_0 - pt.immune_cells);
            pt.muscle_cells += dt_bio * pt.muscle_cells
                * (beta1_smc * eps_str + beta2_smc * eps_elas + beta3_smc * eps_imm);

            pt.fibroblast += dt_bio * ((r_f1 + r_f2 * pt.active_tgf_beta) * pt.fibroblast
                                        - r_f3 * pt.fibroblast);
            pt.procollagen += dt_bio * ((r_p1 + r_p2 * pt.active_tgf_beta) * pt.fibroblast
                                        - r_p3 * pt.procollagen);
            pt.zymogen += dt_bio * ((r_z1 / (1.0 + r_z2 * pt.active_tgf_beta)) * pt.fibroblast
                                        - r_z3 * pt.zymogen);
            pt.collagenase_fib += dt_bio * (r_ca1 * pt.zymogen
                                        - (r_ca2 + r_ca3 * pt.timp) * pt.collagenase_fib);
            pt.timp += dt_bio * ((r_i1 + r_i2 * pt.active_tgf_beta) * pt.fibroblast
                                        - (r_i3 + r_i4 * pt.collagenase_fib) * pt.timp);
            pt.collagen_ad += dt_bio * (r_c1 * pt.procollagen
                                        - r_c2 * pt.collagenase_fib * pt.collagen_ad);

            double term1 = (r_betal1 * pt.active_tgf_beta + r_betal2 * f_lam)
                         / (1.0 + r_betal3 * pt.collagen_ad);
            double term3 = r_betal4 + r_betal5 * f_lam * pt.fibroblast;
            pt.latent_tgf_beta += dt_bio * (term1 * pt.fibroblast * tgf_beta_level
                                           - term3 * pt.latent_tgf_beta);
            pt.active_tgf_beta += dt_bio
                * ((r_beta1 + r_beta2 * f_lam * pt.fibroblast) * pt.latent_tgf_beta
                   - r_beta3 * pt.active_tgf_beta);

            double lc_max  = lam / pt.lambda_rec_min;
            double lc_min  = lam / pt.lambda_rec_max;
            double lc_mode = lam / pt.lambda_rec_mode;

            int idx = std::min(pt.hist_count, FEAneurysmPoint::MAX_HIST - 1);
            pt.lc_max_history[idx] = lc_max;
            pt.hist_count++;

            double lam_att_max  = CalcLamAttMax(pt, dt_bio);
            double lam_att_min  = lam_att_max - width_att_dist;
            double lam_att_mode = lam_att_min + skew_att_dist * (lam_att_min - lam_att_max);

            double col_safe = std::max(0.01, pt.collagen_ad);
            double fib_safe = std::max(0.01, pt.fibroblast);
            double cas_safe = std::max(0.0,  pt.collagenase_fib);
            double alpha    = alpha_init * (fib_safe / col_safe)
                            * sqrt(col_safe * cas_safe);

            alpha = std::min(alpha, 5.0);

            auto safe_div = [](double num, double den) -> double {
                return (std::fabs(den) > 1e-12) ? num / den : 0.0;
            };

            pt.lambda_rec_min  += dt_bio * alpha
                * safe_div(lc_max - lam_att_max, lam_att_max);
            pt.lambda_rec_max  += dt_bio * alpha
                * safe_div(lc_min - lam_att_min, lam_att_min);
            pt.lambda_rec_mode += dt_bio * alpha
                * safe_div(lc_mode - lam_att_mode, lam_att_mode);


            pt.lambda_rec_min  = std::max(0.5, pt.lambda_rec_min);
            pt.lambda_rec_max  = std::max(pt.lambda_rec_min + 0.01, pt.lambda_rec_max);
            pt.lambda_rec_mode = std::max(pt.lambda_rec_min,
                                 std::min(pt.lambda_rec_max, pt.lambda_rec_mode));

            pt.elastin_me = std::max(0.01, pt.elastin_me);
            pt.collagen_me  = std::max(0.01, pt.collagen_me);
            pt.collagen_ad = std::max(0.01, pt.collagen_ad);
            pt.muscle_cells = std::max(0.01, pt.muscle_cells);
            pt.fibroblast = std::max(0.01, pt.fibroblast);
            pt.collagenase_fib = std::max(0.0,  pt.collagenase_fib);
            pt.elastases_imm = std::max(0.0,  pt.elastases_imm);
            pt.collagenases_imm = std::max(0.0,  pt.collagenases_imm);
            pt.zymogen = std::max(0.0,  pt.zymogen);
            pt.timp = std::max(0.0,  pt.timp);
            pt.latent_tgf_beta = std::max(0.0,  pt.latent_tgf_beta);
            pt.active_tgf_beta = std::max(0.0,  pt.active_tgf_beta);
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

        if (t_febio <= febio_disease_start) {
            return c_pass_raw + c_act_raw;
        }

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
    ADD_PROPERTY(m_pActive, "active_solid", FEProperty::Optional);
    ADD_PARAMETER(layer,               "layer");
    ADD_PARAMETER(genotype,            "genotype");
    ADD_PARAMETER(polygenic_score,     "polygenic_score");
    ADD_PARAMETER(febio_disease_start, "febio_disease_start");
    ADD_PARAMETER(years_per_febio_sec, "years_per_febio_sec");
    ADD_PARAMETER(patch_x,             "patch_x");
    ADD_PARAMETER(patch_y,             "patch_y");
    ADD_PARAMETER(patch_z,             "patch_z");
    ADD_PARAMETER(spread_sigma,        "spread_sigma");
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