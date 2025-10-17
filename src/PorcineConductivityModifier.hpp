#ifndef PORCINE_CONDUCTIVITY_MODIFIER_HPP_
#define PORCINE_CONDUCTIVITY_MODIFIER_HPP_

/*
	Element-wise conductivity scaling for ischemic regions in a porcine bidomain setup.
 
 	Behavior:
 	- Constructor takes an ischemia flag, a scaling factor (factor_h), a spherical region
     defined by center (xx,yy,zz) and radius, a "transmural" domain gate (unsigned),
     and a reference to the BidomainProblem<3>.
   	- rCalculateModifiedConductivityTensor(elementIndex, rOriginalConductivity, domainIndex)
     returns a reference to an internal tensor (mTensor) that equals the original tensor
     OR a scaled version when conditions match.
 
  	Notes:
   	- Domain gating: if _transmural==0 -> apply anywhere; if >0 -> apply ONLY when
     domainIndex == _transmural
   	- If the element/node location is inside the ischemic sphere and domain gating passes,
     we multiply the whole tensor by factor_h
   	- The returned tensor preserves symmetry if the original is symmetric.
 */

#include "AbstractConductivityModifier.hpp"
#include "BidomainProblem.hpp"
#include "UblasIncludes.hpp"
#include <cmath>
#include <algorithm>

class PorcineConductivityModifier : public AbstractConductivityModifier<3,3>
{
private:
    // Returned storage
    c_matrix<double,3,3> mTensor;

    // Config
    bool   mIschaemia = false;
    double mFactor    = 1.0;   // multiplicative factor for conductivity
    double mRadius    = 0.0;   // ischemic sphere radius (same units as mesh coordinates)
    double mCx        = 0.0;   // center x
    double mCy        = 0.0;   // center y
    double mCz        = 0.0;   // center z
    unsigned mTransmuralGate = 0; // 0: apply anywhere; >0: apply only if domainIndex==this

    BidomainProblem<3>& mrProblem;

public:
    /*
     * ischaemia   Enable/disable ischemic scaling.
     * factor_h    Multiplicative conductivity factor inside region (e.g., 0.5 to halve).
     * radius      Sphere radius.
     * xx,yy,zz    Sphere center coordinates.
     * transmural  0 => no domain gating (apply anywhere).
     *            >0 => apply ONLY if domainIndex passed to the API equals this value.
     * rProblem    Reference to the bidomain problem (access to mesh/locations).
     */
    PorcineConductivityModifier(bool ischaemia,
                                /*bool phase1a,*/
                                /*bool phase1b,*/
                                double factor_h,
                                double radius,
                                double xx,
                                double yy,
                                double zz,
                                unsigned transmural,
                                BidomainProblem<3>& rProblem)
        : AbstractConductivityModifier<3,3>(),
          mIschaemia(ischaemia),
          mFactor(factor_h),
          mRadius(std::max(0.0, radius)),
          mCx(xx), mCy(yy), mCz(zz),
          mTransmuralGate(transmural),
          mrProblem(rProblem)
    {

    }

    // Setters (optional)
    void SetIschaemia(bool on)             { mIschaemia = on; }
    void SetFactor(double f)               { mFactor = f; }
    void SetSphere(double cx, double cy, double cz, double r)
    {
        mCx = cx; mCy = cy; mCz = cz; mRadius = std::max(0.0, r);
    }
    void SetTransmuralGate(unsigned gate)  { mTransmuralGate = gate; }

    /*
     * Return a possibly modified conductivity tensor for the given element index.
     *
     * elementIndex         	Index of the element
     * rOriginalConductivity 	Original conductivity tensor 
     * domainIndex          	Domain label
     */
    c_matrix<double,3,3>& rCalculateModifiedConductivityTensor(unsigned elementIndex,
                                                               const c_matrix<double,3,3>& rOriginalConductivity,
                                                               unsigned domainIndex)
    {
        // Start from original
        mTensor = rOriginalConductivity;

        // Quick exits if nothing to do
        if (!mIschaemia)            { return mTensor; }
        if (mRadius <= 0.0)         { return mTensor; }
        if (mFactor == 1.0)         { return mTensor; }

        // Domain gate: if mTransmuralGate > 0, apply only within that domain
        if (mTransmuralGate > 0 && domainIndex != mTransmuralGate)
        {
            return mTensor;
        }

        // Get spatial location.
        const double x = mrProblem.rGetMesh().GetNode(elementIndex)->rGetLocation()[0];
        const double y = mrProblem.rGetMesh().GetNode(elementIndex)->rGetLocation()[1];
        const double z = mrProblem.rGetMesh().GetNode(elementIndex)->rGetLocation()[2];

        const double dx = x - mCx;
        const double dy = y - mCy;
        const double dz = z - mCz;
        const double r2 = dx*dx + dy*dy + dz*dz;

        if (r2 <= mRadius*mRadius)
        {
            // Inside ischemic sphere: scale full tensor.
            mTensor = mTensor * mFactor;
        }

        return mTensor;
    }
};

#endif // PORCINE_CONDUCTIVITY_MODIFIER_HPP_
