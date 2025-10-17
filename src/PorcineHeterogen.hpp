#ifndef PORCINEHETEROGEN_HPP_
#define PORCINEHETEROGEN_HPP_

/*
  Region-wise conductance scaling to impose porcine ventricular heterogeneities.
 
  Usage:
  PorcineHeterogen hetero(&heart_geom_info, 
  							apexBase=true, transmural=true, interventricular=true)
    
   hetero.SetApexBaseAnchors(apex_z, base_z); // <-- reuse your previous numbers here

 */

#include <cmath>
#include <algorithm>

// Forward declarations
template<unsigned ELEMENT_DIM> class Node;
template<unsigned SPACE_DIM> class HeartGeometryInformation;

class PorcineHeterogen
{
public:
    PorcineHeterogen(const HeartGeometryInformation<3>* pHeartGeomInfo,
                     bool enableApexBase,
                     bool enableTransmural,
                     bool enableInterventricular)
        : mpHeartGeomInfo(pHeartGeomInfo),
          mApexBase(enableApexBase),
          mTransmural(enableTransmural),
          mInterventricular(enableInterventricular)
    {
        mNormaliser = (mApexBase ? 1 : 0) + (mTransmural ? 1 : 0) + (mInterventricular ? 1 : 0);
        if (mNormaliser == 0) { mNormaliser = 1; }
    }

	// Allow wiring geometry after construction.
    void SetGeometry(const HeartGeometryInformation<3>* pHeartGeomInfo)
    {
        mpHeartGeomInfo = pHeartGeomInfo;
    }

	// Config
    void SetApexBaseAnchors(double apex_z, double base_z)
    {
        mHasApexBaseAnchors = true;
        mApexZ = apex_z;
        mBaseZ = base_z;
        // Precompute mid & half-span for speed; guard zero span.
        mAbMid      = 0.5*(mApexZ + mBaseZ);
        mAbHalfSpan = std::max(std::abs(mBaseZ - mApexZ) * 0.5, 1e-12);
    }

    /// Set heterogeneity multipliers
    void SetMultipliers(double apexBaseMultiplier,
                        double transmuralMultiplier,
                        double interventricularMultiplier)
    {
        mApexBaseMultiplier         = ClampMultiplier(apexBaseMultiplier);
        mTransmuralMultiplier       = ClampMultiplier(transmuralMultiplier);
        mInterventricularMultiplier = ClampMultiplier(interventricularMultiplier);
    }

    void SetOutputClamp(double min_scale, double max_scale)
    {
        mOutMin = std::min(min_scale, max_scale);
        mOutMax = std::max(min_scale, max_scale);
    }

    // Main

    /*
     * Compute per-node scaling (unitless), averaged over enabled axes and clamped.
     */
    double GetScalingFactor(const Node<3>& rNode) const

	if (!mpHeartGeomInfo)
    {
        // Geometry not wired yet; return neutral scaling
        return 1.0;
    }
	
    {
        double scaling_sum = 0.0;

        // --- Apex–base ---
        if (mApexBase)
        {
            double s_ab = 1.0;
            if (mHasApexBaseAnchors)
            {
                const double z = rNode.rGetLocation()[2];
                // Signed projection in [-1,+1]: + toward base_z, - toward apex_z.
                const double proj = ClampUnit((z - mAbMid) / mAbHalfSpan);
                s_ab = BipolarScale(proj, mApexBaseMultiplier);
            }
            // else neutral (1.0) if not set
            scaling_sum += s_ab;
        }

        // --- Transmural ---
        if (mTransmural)
        {
            const auto idx = rNode.GetIndex();
            const double d_endo = mpHeartGeomInfo->rGetDistanceMapEndocardium()[idx];
            const double d_epi  = mpHeartGeomInfo->rGetDistanceMapEpicardium()[idx];
            const double proj_tm = SignedRatio(d_epi, d_endo); // in [-1,+1]
            const double s_tm = BipolarScale(proj_tm, mTransmuralMultiplier);
            scaling_sum += s_tm;
            //scaling_sum += 1.0;
        }

        // --- Interventricular ---
        if (mInterventricular)
        {
            const unsigned idx = rNode.GetIndex();
            const double lv_d  = mpHeartGeomInfo->rGetDistanceMapLeftVentricle()[idx];
            const double rv_d  = mpHeartGeomInfo->rGetDistanceMapRightVentricle()[idx];

            double s_iv = 1.0;
            if (lv_d < rv_d)       { s_iv = mInterventricularMultiplier; }
            else if (rv_d < lv_d)  { s_iv = 1.0 / mInterventricularMultiplier; }
            scaling_sum += s_iv;
        }

        const double avg = scaling_sum / static_cast<double>(mNormaliser);
        return std::max(mOutMin, std::min(mOutMax, avg));
    }

private:
    // helpers
    static double ClampMultiplier(double v)
    {
        const double kMin = 1.01; // preserve at least 1% effect if user sets too low
        const double kMax = 1.50; // cap at +50% 
        if (v < kMin) return kMin;
        if (v > kMax) return kMax;
        return v;
    }

    static double ClampUnit(double v)
    {
        if (v < -1.0) return -1.0;
        if (v >  1.0) return  1.0;
        return v;
    }

    // Map signed projection [-1,+1] to multiplier (>1) or its reciprocal (<1).
    static double BipolarScale(double signed_projection, double multiplier)
    {
        if (signed_projection > 0.0) return multiplier;        // epi
        if (signed_projection < 0.0) return 1.0 / multiplier;  // endo
        return 1.0;
    }

    // Convert two positive distances to a signed ratio in [-1,+1].
    static double SignedRatio(double pos_dist_plus, double pos_dist_minus)
    {
        const double eps = 1e-12;
        const double num = pos_dist_plus - pos_dist_minus;
        const double den = pos_dist_plus + pos_dist_minus + eps;
        return ClampUnit(num / den);
    }

private:
    const HeartGeometryInformation<3>* mpHeartGeomInfo;

    bool mApexBase;
    bool mTransmural;
    bool mInterventricular;
    int  mNormaliser;

    // Apex-base anchors
    bool   mHasApexBaseAnchors = false;
    double mApexZ = 0.0;
    double mBaseZ = 0.0;
    double mAbMid = 0.0;
    double mAbHalfSpan = 1.0;

    // Tunable multipliers
    double mApexBaseMultiplier         = 1.20;
    double mTransmuralMultiplier       = 1.15;
    double mInterventricularMultiplier = 1.10;

    // Final clamp
    double mOutMin = 0.70;
    double mOutMax = 1.30;
};

#endif // PORCINEHETEROGEN_HPP_
