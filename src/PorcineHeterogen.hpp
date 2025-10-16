#ifndef PORCINEHETEROGEN_HPP_
#define PORCINEHETEROGEN_HPP_

#include "HeartGeometryInformation.hpp"

class PorcineHeterogen
{
private:
    HeartGeometryInformation<3>* mpHeartGeomInfo;
    bool mApexBase;
    bool mTransmural;
    bool mInterventr;
    int mNormaliser;

public:
    
    PorcineHeterogen(HeartGeometryInformation<3>* pHeartGeomInfo, bool ApexBase, bool Transmural, bool Interventr) :
		mpHeartGeomInfo(pHeartGeomInfo),
		mApexBase(ApexBase),
		mTransmural(Transmural),
		mInterventr(Interventr)
    {
        mNormaliser = (int)mApexBase+(int)mTransmural+(int)mInterventr;
    }
    
	double GetScalingFactor(const Node<3>& rNode)
	{
		if (mNormaliser==0) return 1.0;

		unsigned node_index = rNode.GetIndex();
		HeartRegionType this_region = mpHeartGeomInfo->GetHeartRegion(node_index);
		double scaling = 0.0;
		double multiplier = 0.2;
		
		if (mApexBase)
		{
			const double z = rNode.GetPoint()[2];
			const double zmin = 47.246; // apex coord
			const double zmax = 54.383; // base coord
			const double width = zmax-zmin;
			const double coord = 2.0*(z-zmin)/width - 1.0; // current coord
			
			scaling += pow(multiplier, coord);
		}
		
		if (mTransmural)
		{
			const double coord = 1.0-2.0*mpHeartGeomInfo->CalculateRelativeWallPosition(node_index);
			scaling += pow(multiplier, coord);
		}
		
		if (mInterventr)
		{
			if ( this_region == mpHeartGeomInfo->RIGHT_VENTRICLE_WALL )
			{
				scaling += 1.0/multiplier;
			}
			else
			{
				scaling += multiplier;
			}
		}
		return scaling/mNormaliser;
	}

};

#endif
