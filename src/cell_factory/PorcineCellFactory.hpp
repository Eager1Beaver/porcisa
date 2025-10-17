#ifndef PORCINE_CELL_FACTORY_HPP_
#define PORCINE_CELL_FACTORY_HPP_

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "UblasIncludes.hpp"

#include "AbstractCardiacCellFactory.hpp"
#include "HeartGeometryInformation.hpp"
#include "AbstractCvodeCell.hpp"
#include "AbstractCardiacCellInterface.hpp"
#include "SimpleStimulus.hpp"
#include "ChasteEllipsoid.hpp"

#include "pig_ventr_ap_epiCvodeOpt.hpp"
#include "pig_ventr_ap_endoCvodeOpt.hpp"
#include "pig_ventr_ap_mCvodeOpt.hpp"

// forward
class PorcineHeterogen;

/*
	Cell factory with setters for layers, ischemia, and multi-region stimulation.
 */
class PorcineCellFactory : public AbstractCardiacCellFactory<3>
{
public:
    // ---- Construction ----
    explicit PorcineCellFactory(
        boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
        PorcineHeterogen* pHetero,
        bool ischaemia = false)
        : AbstractCardiacCellFactory<3>(pSolver),
          mpPorcineHeterogen(pHetero),
          mIschaemia(ischaemia)
    {
        // Defaults mirroring typical porcine transmural composition
        SetLayerFractionsFromPercentages(30.0, 40.0, 30.0); // epi/m/endo in %
        // Default transmural bands
        SetTransmuralBands(/*band0_low=*/0.1, /*band0_high=*/0.5,
                           /*band1_low=*/0.1, /*band1_high=*/1.0, // (>0.1)
                           /*band2_low=*/0.1, /*band2_high=*/0.3);
        // Default ischemia geometry disabled until set
        SetIschemiaGeometry(/*cx=*/0, /*cy=*/0, /*cz=*/0, /*icz_r=*/0, /*bz_r=*/0);
        // Default ischemia multipliers
        SetIschemiaMultipliersCore(1,1,1,1,1,1,1);
        SetIschemiaMultipliersBorder(1,1,1,1,1,1,1);
        SetBorderLinearTaper(true);
        // Stimulus defaults
        mStimAmplitude = -45000.0; // uA/cm^3
        mStimDurationMs = 2.0;
    }

    explicit PorcineCellFactory(
        boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
        PorcineHeterogen* pHetero,
        const std::vector<double>& stimulus_starts_ms,
        bool ischaemia = false)
        : PorcineCellFactory(pSolver, pHetero, ischaemia)
    {
        mStimulusStartsMs = stimulus_starts_ms;
    }

    // ---- Layering controls ----
    void SetLayerFractions(double epi_frac, double m_frac, double endo_frac)
    {
        const double sum = epi_frac + m_frac + endo_frac;
        if (sum <= 0) throw std::runtime_error("Layer fractions must be positive");
        mLayerEpi  = epi_frac  / sum;
        mLayerMid  = m_frac    / sum;
        mLayerEndo = endo_frac / sum;
    }
    void SetLayerFractionsFromPercentages(double epi_pct, double m_pct, double endo_pct)
    { SetLayerFractions(epi_pct/100.0, m_pct/100.0, endo_pct/100.0); }

    // ---- Ischemia geometry: ICZ (core) radius and BZ (border) radius > ICZ ----
    void SetIschemiaGeometry(double cx, double cy, double cz, double icz_radius, double bz_radius)
    {
        mCx = cx; mCy = cy; mCz = cz;
        mIczR = std::max(0.0, icz_radius);
        mBzR  = std::max(mIczR, bz_radius); // ensure BZ >= ICZ
    }

    // ---- Transmural gating bands ----
    // Band 0: default gate used when mTransmuralMode==0  -> (low0, high0)
    // Band 1: when mTransmuralMode==1                    -> (low1, high1)
    // Band 2: when mTransmuralMode==2                    -> (low2, high2)
    void SetTransmuralBands(double low0, double high0,
                            double low1, double high1,
                            double low2, double high2)
    {
        mBand0 = {low0, high0};
        mBand1 = {low1, high1};
        mBand2 = {low2, high2};
    }
    void SetTransmuralMode(unsigned mode) { mTransmuralMode = mode; }

    // ---- Ischemia multipliers (unitless); core vs border; 1.0 means no change ----
    // Order: gK_o (extracellular K), gNa, gCaL, gKATP, gNCX, PNaK, gNaL
    void SetIschemiaMultipliersCore(double Ko, double gNa, double gCaL, double gKATP, double gNCX, double PNaK, double gNaL)
    { mCore = {Ko,gNa,gCaL,gKATP,gNCX,PNaK,gNaL}; }
    void SetIschemiaMultipliersBorder(double Ko, double gNa, double gCaL, double gKATP, double gNCX, double PNaK, double gNaL)
    { mBorder = {Ko,gNa,gCaL,gKATP,gNCX,PNaK,gNaL}; }
    void SetBorderLinearTaper(bool on) { mBorderLinearTaper = on; }

    // ---- Stimulus controls ----
    // Supply ellipsoid regions and starts (ms). Amplitudes/duration are global setters.
    void SetStimulusRegions(const std::vector<c_vector<double,3>>& centers,
                            const std::vector<double>& radii,
                            const std::vector<double>& starts_ms)
    {
        if (centers.size() != radii.size() || centers.size() != starts_ms.size())
            throw std::runtime_error("Stimulus vectors must have equal length");
        mStimCenters = centers;
        mStimRadii   = radii;
        mStimulusStartsMs = starts_ms;
    }

	// ---- Define A..H by centers/radii ----
    void SetEightRootStimuli(const std::vector<c_vector<double,3>>& centers,
                             const std::vector<double>& radii,
                             double start0_ms,
                             double step_msA=3.5, double step_msB=6.0, double step_msC=16.0,
                             double step_msD=8.0,  double step_msE=18.0, double step_msF=20.0, double step_msG=22.0)
    {
        if (centers.size() != 8 || radii.size() != 8)
            throw std::runtime_error("SetEightRootStimuli expects 8 centers and 8 radii");
        std::vector<double> starts = {
            start0_ms,
            start0_ms + step_msA,
            start0_ms + step_msB,
            start0_ms + step_msC,
            start0_ms + step_msD,
            start0_ms + step_msE,
            start0_ms + step_msF,
            start0_ms + step_msG
        };
        SetStimulusRegions(centers, radii, starts);
    }
	
	/// Load user's eight spherical root-stimuli (radius 0.25) with representative start times.
	void UseUserRootStimulusPreset()
    {
        std::vector<c_vector<double,3>> centers(8);
        std::vector<double> radii(8, 0.25);
        std::vector<double> starts_ms = {
            5.0,   // A: basal septum (5–17 ms)
            7.0,   // B: central septum (7–10 ms)
            11.0,  // C: anterior apex (11–15–18 ms)
            23.0,  // D: posterior central (~23 ms)
            15.0,  // E: posterior apex (15–20 ms)
            15.0,  // F: anterior central (~15 ms)
            22.0,  // G: LV free wall central (~22 ms)
            11.0   // H: RV posterior apex (11–15 ms)
        };

        // A
        centers[0][0] = 21.8 - 0.264; centers[0][1] = 6.9 + 0.793; centers[0][2] = 54.5 - 0.83  - 0.3;
        // B
        centers[1][0] = 21.8 - 0.21;  centers[1][1] = 6.9 - 0.59;  centers[1][2] = 54.5 - 2.46  - 0.3;
        // C
        centers[2][0] = 21.8 + 1.36;  centers[2][1] = 6.9 - 2.55;  centers[2][2] = 54.5 - 4.90  - 0.2;
        // D
        centers[3][0] = 21.8 + 1.40;  centers[3][1] = 6.9 + 2.07;  centers[3][2] = 54.5 - 3.35  - 0.3;
        // E
        centers[4][0] = 21.8 + 1.84;  centers[4][1] = 6.9 + 1.09;  centers[4][2] = 54.5 - 5.25  - 0.2;
        // F
        centers[5][0] = 21.8 + 2.25;  centers[5][1] = 6.9 - 2.16;  centers[5][2] = 54.5 - 2.60  - 0.3;
        // G
        centers[6][0] = 21.8 + 3.60;  centers[6][1] = 6.9 + 0.66;  centers[6][2] = 54.5 - 1.58  - 0.3;
        // H
        centers[7][0] = 21.8 - 1.70;  centers[7][1] = 6.9 - 2.38;  centers[7][2] = 54.5 - 2.626 - 0.3;

        SetStimulusRegions(centers, radii, starts_ms);
    }

    void SetStimulusAmplitude(double a) { mStimAmplitude = a; }
    void SetStimulusDurationMs(double d_ms) { mStimDurationMs = d_ms; }

    // Store starts only
    void SetStimulusStarts(const std::vector<double>& starts_ms) { mStimulusStartsMs = starts_ms; }    

    // ---- Main factory ----
    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<3>* pNode) override
    {
        AbstractCvodeCell* p_cell = nullptr;

        const unsigned node_index = pNode->GetIndex();
        HeartGeometryInformation<3>* p_geom = this->GetHeartGeometryInformation();

        // --- Transmural relative position: endo(0) -> epi(1)
        const double lv_d = p_geom->rGetDistanceMapLeftVentricle()[node_index];
        const double rv_d = p_geom->rGetDistanceMapRightVentricle()[node_index];
        const double d_endo = std::min(lv_d, rv_d);
        const double d_epi  = p_geom->rGetDistanceMapEpicardium()[node_index];
        const double rel_pos = d_endo/(d_endo + d_epi + 1e-12); // [0,1]

        // --- Cell type by layer (epi / mid / endo)
        if (rel_pos > (1.0 - mLayerEpi))
        {
            p_cell = new pig_ventr_ap_epiCvodeOpt(this->mpSolver, this->mpZeroStimulus);
        }
        else if (rel_pos > (1.0 - mLayerEpi - mLayerMid))
        {
            p_cell = new pig_ventr_ap_mCvodeOpt(this->mpSolver, this->mpZeroStimulus);
        }
        else
        {
            p_cell = new pig_ventr_ap_endoCvodeOpt(this->mpSolver, this->mpZeroStimulus);
        }

        // Ischemia effects
        if (mIschaemia && (mIczR > 0.0 || mBzR > 0.0))
        {
            const double x = pNode->rGetLocation()[0];
            const double y = pNode->rGetLocation()[1];
            const double z = pNode->rGetLocation()[2];
            const double dx = x - mCx, dy = y - mCy, dz = z - mCz;
            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            const bool in_core  = (mIczR > 0.0) && (r <= mIczR);
            const bool in_bz    = (mBzR  > mIczR) && (r > mIczR) && (r <= mBzR);

            // Transmural gate
            if (PassesTransmuralGate(rel_pos))
            {
                if (in_core)
                {
                    ApplyIschemiaMultipliers(*p_cell, mCore);
                }
                else if (in_bz)
                {
                    if (mBorderLinearTaper)
                    {
                        const double t = (r - mIczR) / std::max(1e-12, (mBzR - mIczR)); // 0..1 across BZ
                        IschParams tapered = Lerp(mCore, mBorder, t);
                        ApplyIschemiaMultipliers(*p_cell, tapered);
                    }
                    else
                    {
                        ApplyIschemiaMultipliers(*p_cell, mBorder);
                    }
                }
            }
        }

        // --- Region-wise stimuli (optional)
        if (!mStimCenters.empty())
        {
            // build ellipsoids and stimuli
            EnsureStimulusObjects();

            // assign first matching region’s stimulus (or you can combine if needed)
            const c_vector<double,3>& xyz = pNode->rGetLocation();
            for (size_t i = 0; i < mEllipsoids.size(); ++i)
            {
                if (mEllipsoids[i].DoesContain(xyz))
                {
                    p_cell->SetStimulusFunction(mSimpleStimuli[i]);
                    break;
                }
            }
        }
        else if (!mStimulusStartsMs.empty())
        {
            // If starts were provided but no regions, use the first start as a global stimulus
            boost::shared_ptr<SimpleStimulus> p_stim(new SimpleStimulus(mStimAmplitude, mStimDurationMs, mStimulusStartsMs.front()));
            p_cell->SetStimulusFunction(p_stim);
            // keep a reference to avoid re-allocating identical stimuli
            if (mSimpleStimuli.empty()) mSimpleStimuli.push_back(p_stim);
        }

        // --- Heterogeneity scaling for IKs
        const double gks = p_cell->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
        const double scale = (mpPorcineHeterogen ? mpPorcineHeterogen->GetScalingFactor(*pNode) : 1.0);
        p_cell->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", scale * gks);

        // Tight tolerances
        p_cell->SetTolerances(1e-4, 1e-6);

        return p_cell;
    }

private:
    // ---- internal structures ----
    struct IschParams { double Ko, gNa, gCaL, gKATP, gNCX, PNaK, gNaL; };

    // ---- helpers ----
    static double Clamp01(double x) { return x < 0 ? 0 : (x > 1 ? 1 : x); }

    static IschParams Lerp(const IschParams& a, const IschParams& b, double t)
    {
        t = Clamp01(t);
        return {
            a.Ko   + (b.Ko   - a.Ko)   * t,
            a.gNa  + (b.gNa  - a.gNa)  * t,
            a.gCaL + (b.gCaL - a.gCaL) * t,
            a.gKATP+ (b.gKATP- a.gKATP)* t,
            a.gNCX + (b.gNCX - a.gNCX) * t,
            a.PNaK + (b.PNaK - a.PNaK) * t,
            a.gNaL + (b.gNaL - a.gNaL) * t
        };
    }

    bool PassesTransmuralGate(double rel_pos) const
    {
        const std::pair<double,double>* band = nullptr;
        if      (mTransmuralMode == 0) band = &mBand0;
        else if (mTransmuralMode == 1) band = &mBand1;
        else                            band = &mBand2;
        return (rel_pos > band->first) && (rel_pos < band->second);
    }

    void ApplyIschemiaMultipliers(AbstractCvodeCell& cell, const IschParams& P) const
    {
        // K_o
        cell.SetParameter("extracellular_potassium_concentration",
                          P.Ko * cell.GetParameter("extracellular_potassium_concentration"));

        // gNa, gCaL, gKATP
        const double gna = cell.GetParameter("membrane_fast_sodium_current_conductance");
        const double gca = cell.GetParameter("membrane_L_type_calcium_current_conductance");
        const double gkatp = cell.GetParameter("membrane_atp_dependent_potassium_current_conductance");
        cell.SetParameter("membrane_fast_sodium_current_conductance", gna * P.gNa);
        cell.SetParameter("membrane_L_type_calcium_current_conductance", gca * P.gCaL);
        cell.SetParameter("membrane_atp_dependent_potassium_current_conductance", gkatp * P.gKATP);

        // NCX
        const double gNCX = cell.GetParameter("membrane_sodium_calcium_exchanger_current_conductance");
        cell.SetParameter("membrane_sodium_calcium_exchanger_current_conductance", gNCX * P.gNCX);

        // NaK
        const double pNaK = cell.GetParameter("membrane_sodium_potassium_pump_current_permeability");
        cell.SetParameter("membrane_sodium_potassium_pump_current_permeability", pNaK * P.PNaK);

        // NaL
        const double gNaL = cell.GetParameter("membrane_late_sodium_current_conductance");
        cell.SetParameter("membrane_late_sodium_current_conductance", gNaL * P.gNaL);
    }

    void EnsureStimulusObjects()
    {
        if (!mEllipsoids.empty()) return;
        mEllipsoids.reserve(mStimCenters.size());
        mSimpleStimuli.reserve(mStimCenters.size());

        for (size_t i = 0; i < mStimCenters.size(); ++i)
        {
            const c_vector<double,3> c = mStimCenters[i];
            const double r = mStimRadii[i];
            c_vector<double,3> radii = zero_vector<double>(3);
            radii[0] = r; radii[1] = r; radii[2] = r;
            mEllipsoids.emplace_back(c, radii);

            mSimpleStimuli.push_back(
                boost::shared_ptr<SimpleStimulus>(
                    new SimpleStimulus(mStimAmplitude, mStimDurationMs, mStimulusStartsMs[i])
                )
            );
        }
    }

private:
    // ---- heterogeneity ----
    PorcineHeterogen* mpPorcineHeterogen = nullptr;

    // ---- layer fractions ----
    double mLayerEpi  = 0.3;
    double mLayerMid  = 0.4;
    double mLayerEndo = 0.3;

    // ---- ischemia geometry ----
    double mCx=0, mCy=0, mCz=0, mIczR=0, mBzR=0;
    bool   mIschaemia = false;

    // transmural gating
    unsigned mTransmuralMode = 0;
    std::pair<double,double> mBand0{0.1,0.5};
    std::pair<double,double> mBand1{0.1,1.0};
    std::pair<double,double> mBand2{0.1,0.3};

    // ischemia multipliers
    IschParams mCore   {1,1,1,1,1,1,1};
    IschParams mBorder {1,1,1,1,1,1,1};
    bool mBorderLinearTaper = true;

    // ---- stimuli ----
    double mStimAmplitude = -45000.0;
    double mStimDurationMs = 2.0;
    std::vector<double> mStimulusStartsMs;
    std::vector<c_vector<double,3>> mStimCenters;
    std::vector<double> mStimRadii;
    std::vector<ChasteEllipsoid<3>> mEllipsoids;
    std::vector<boost::shared_ptr<SimpleStimulus>> mSimpleStimuli;
};

#endif // PORCINE_CELL_FACTORY_HPP_
