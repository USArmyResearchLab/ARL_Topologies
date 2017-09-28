/*
 * ARL_Topologies - An extensible topology optimization program
 * 
 * Written in 2017 by Raymond A. Wildman <raymond.a.wildman.civ@mail.mil>
 * This project constitutes a work of the United States Government and is not 
 * subject to domestic copyright protection under 17 USC Sec. 105.
 * Release authorized by the US Army Research Laboratory
 * 
 * To the extent possible under law, the author(s) have dedicated all copyright 
 * and related and neighboring rights to this software to the public domain 
 * worldwide. This software is distributed without any warranty.
 * 
 * You should have received a copy of the CC0 Public Domain Dedication along 
 * with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 
 * 
 */

// Questions?
// Contact: Raymond Wildman, raymond.a.wildman.civ@mail.mil

#include "torfactory.h"
#include "topoptrep.h"
#include "pixelrep.h"
#include "volmesh2d.h"
#include "volmesh3d.h"
#include "voxelrep.h"
#include "csgtree.h"
#include "inputloaderrep.h"
#include "helper.h"

namespace Topologies{
namespace TopOptRepFactory
{
	std::unique_ptr<TopOptRep> createTopOptRep(TORType inTORT, const InputLoader::OptNodeInfo& oni)
	{
		// Construct RNI from the ONI
		InputLoader::RepNodeInfo rni(inTORT, oni);
		return createTopOptRep(rni);
	}

	std::unique_ptr<TopOptRep> createTopOptRep(const InputLoader::RepNodeInfo& rni)
	{
		std::unique_ptr<TopOptRep> upTOR;
		std::vector<double> penalParams;
		std::vector<double> projParams;
		if(rni.getType() == tortPixel)
		{
			InputLoader::TORGenericVolume torParser(rni.getTypeName());
			torParser.parseNode(rni);
			penalParams = {torParser.getPenalPower(), torParser.getMinDensity()};
			upTOR = std::unique_ptr<TopOptRep>(new PixelRep<>(tortPixel, torParser, penalParams, projParams));
		}
		else if(rni.getType() == tortMesh2D)
		{
			InputLoader::TORGenericMesh torParser(rni.getTypeName());
			torParser.parseNode(rni);
			penalParams = {torParser.getPenalPower(), torParser.getMinDensity()};
			upTOR = std::unique_ptr<TopOptRep>(new VolMesh2D<>(tortMesh2D, torParser, penalParams, projParams));
		}
		else if(rni.getType() == tortHeaviside2D)
		{
			InputLoader::TORGenericVolume torParser(rni.getTypeName());
			torParser.parseNode(rni);
			penalParams = {torParser.getPenalPower(), torParser.getMinDensity()};
			projParams = {torParser.getThreshold(), torParser.getBetaHeavi()};
			upTOR = std::unique_ptr<TopOptRep>(new PixelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>(
					tortHeaviside2D, torParser, penalParams, projParams));
		}
		else if(rni.getType() == tortHeavisideMesh2D)
		{
			InputLoader::TORGenericMesh torParser(rni.getTypeName());
			torParser.parseNode(rni);
			penalParams = {torParser.getPenalPower(), torParser.getMinDensity()};
			projParams = {torParser.getThreshold(), torParser.getBetaHeavi()};
			upTOR = std::unique_ptr<TopOptRep>(new VolMesh2D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>(
					tortHeavisideMesh2D, torParser, penalParams, projParams));
		}
#ifdef USE_EXP
		else if(rni.getType() == tortLowRankPixel)
		{
			InputLoader::TORGenericVolume torParser(rni.getTypeName());
			torParser.parseNode(rni);
			upTOR = std::unique_ptr<TopOptRep>(new LowRankPixelRep(torParser));
		}
#endif
		else if(rni.getType() == tortMesh3D)
		{
			InputLoader::TORGenericMesh torParser(rni.getTypeName());
			torParser.parseNode(rni);
			penalParams = {torParser.getPenalPower(), torParser.getMinDensity()};
			upTOR = std::unique_ptr<TopOptRep>(new VolMesh3D<>(rni.getType(), torParser, penalParams, projParams));
		}
		else if(rni.getType() == tortVoxel)
		{
			InputLoader::TORGenericVolume torParser(rni.getTypeName());
			torParser.parseNode(rni);
			penalParams = {torParser.getPenalPower(), torParser.getMinDensity()};
			upTOR = std::unique_ptr<TopOptRep>(new VoxelRep<>(rni.getType(), torParser, penalParams, projParams));
		}
		else if(rni.getType() == tortHeavisideMesh3D)
		{
			InputLoader::TORGenericMesh torParser(rni.getTypeName());
			torParser.parseNode(rni);
			penalParams = {torParser.getPenalPower(), torParser.getMinDensity()};
			projParams = {torParser.getThreshold(), torParser.getBetaHeavi()};
			upTOR = std::unique_ptr<TopOptRep>(new VolMesh3D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>(rni.getType(),
					torParser, penalParams, projParams));
		}
		else if(rni.getType() == tortHeaviside3D)
		{
			InputLoader::TORGenericVolume torParser(rni.getTypeName());
			torParser.parseNode(rni);
			penalParams = {torParser.getPenalPower(), torParser.getMinDensity()};
			projParams = {torParser.getThreshold(), torParser.getBetaHeavi()};
			upTOR = std::unique_ptr<TopOptRep>(new VoxelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>(rni.getType(),
					torParser, penalParams, projParams));
		}
		else if(rni.getType() == tortCSG2D)
		{
			InputLoader::TORCSGTree torParser(rni.getTypeName());
			torParser.parseNode(rni);
			upTOR = std::unique_ptr<TopOptRep>(new CSGTreeRep(torParser));
		}
		else
		{
			std::cout << "Unknown topology representation, aborting" << std::endl;
			abort();
		}
		return upTOR;
	}

	std::unique_ptr<TopOptRep> createTopOptRep(TORType inTORT, const std::vector<std::vector<double>>& realRep,
			const std::vector<std::vector<int>>& discreteRep)
	{
		std::unique_ptr<TopOptRep> myTOR;
		if(inTORT == tortPixel)
			myTOR = std::unique_ptr<TopOptRep>(new PixelRep<>(inTORT, discreteRep, realRep));
		else if(inTORT == tortMesh2D)
			myTOR = std::unique_ptr<TopOptRep>(new VolMesh2D<>(inTORT, discreteRep, realRep));
		else if(inTORT == tortHeaviside2D)
			myTOR = std::unique_ptr<TopOptRep>(new PixelRep<HelperNS::powPenalMin,HelperNS::thresholdHeaviside>(
					inTORT, discreteRep, realRep));
		else if(inTORT == tortHeavisideMesh2D)
			myTOR = std::unique_ptr<TopOptRep>(new VolMesh2D<HelperNS::powPenalMin,HelperNS::thresholdHeaviside>(
				inTORT, discreteRep, realRep));
		else if(inTORT == tortCSG2D)
			myTOR = std::unique_ptr<TopOptRep>(new CSGTreeRep(discreteRep, realRep));
#ifdef USE_EXP
		else if(inTORT == tortLinearSpline)
			myTOR = std::unique_ptr<TopOptRep>(new LinearSpline(discreteRep, realRep));
		else if(inTORT == tortLowRankPixel)
			myTOR = std::unique_ptr<TopOptRep>(new LowRankPixelRep(discreteRep, realRep));
#endif
		else if(inTORT == tortVoxel)
			myTOR = std::unique_ptr<TopOptRep>(new VoxelRep<>(inTORT, discreteRep, realRep));
		else if(inTORT == tortMesh3D)
			myTOR = std::unique_ptr<TopOptRep>(new VolMesh3D<>(inTORT, discreteRep, realRep));
		else if(inTORT == tortHeaviside3D)
      myTOR = std::unique_ptr<TopOptRep>(new VoxelRep<HelperNS::powPenalMin,HelperNS::thresholdHeaviside>(
          inTORT, discreteRep, realRep));
    else if(inTORT == tortHeavisideMesh3D)
      myTOR = std::unique_ptr<TopOptRep>(new VolMesh3D<HelperNS::powPenalMin,HelperNS::thresholdHeaviside>(
        inTORT, discreteRep, realRep));
		else
		{
			std::cout << "Error, attempting to convert to unknown TOR" << std::endl;
			abort();
		}
		return myTOR;
	}
}
}
