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

#include "volmesh2d.h"
#include "geometrytranslation.h"
#include "helper.h"
#include "filter2d.h"
#include "tomesh.h"
#include "exotxtmeshloader.h"
#include "gmshtxtloader.h"
//#include <random>
#include <algorithm>
#include <unordered_set>
#include "outputwriter.h"
#include "tomeshprocessing.h"

namespace Topologies{
template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh2D<PenaltyFunc, ProjectionFunc>::VolMesh2D(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, 
	const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
	VolMesh<PenaltyFunc, ProjectionFunc>(inTORT, inputParams, inPenalParams, inProjParams)
{
	// Mesh initialization is handled by PixelRep
}


template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh2D<PenaltyFunc, ProjectionFunc>::VolMesh2D(TORType inTORT, const InputLoader::TORGenericMesh& inputParams, 
	const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
	VolMesh<PenaltyFunc, ProjectionFunc>(inTORT, inputParams, inPenalParams, inProjParams)
{
	finishSetup(inputParams.getPolyVec());
}

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh2D<PenaltyFunc, ProjectionFunc>::VolMesh2D(TORType inTORT, const std::vector<std::vector<int> >& discreteParams, 
	const std::vector<std::vector<double> >& realParams) :
	VolMesh<PenaltyFunc, ProjectionFunc>(inTORT)
{
	bool error = realParams.size() < 4;
	std::vector<std::vector<Point_2_base> > ptVecs;
	if(!error)
	{
		error = realParams[0].size() != 4;
		if(!error)
		{
			VM::threshold = realParams[0][0];
			VM::filtRad = realParams[0][1];
			VM::meshParams.triMeshEdgeAngle = realParams[0][2];
			VM::meshParams.triMeshEdgeSize = realParams[0][3];
		}
		// Function parameters
		VM::penalParams = realParams[1];
		VM::projParams = realParams[2];
		// Fixed block parameters
		VM::fixedBlockVec.resize(realParams[3].size());
		for(std::size_t k = 0; k < realParams[3].size(); ++k)
			VM::fixedBlockVec[k].second = realParams[3][k];
		// Boundary points
		std::size_t nPVecs = 4;
		ptVecs.resize(realParams.size() - nPVecs);
		for(std::size_t k = 0; k < ptVecs.size() && !error; ++k)
		{
			error |= (realParams[k + nPVecs].size() % 2) != 0;
			if(!error)
			{
				std::vector<Point_2_base> curPtVec(realParams[k + nPVecs].size()/2);
				for(std::size_t kp = 0; kp < curPtVec.size(); ++kp)
					curPtVec[kp] = Point_2_base(realParams[k + nPVecs][2*kp], realParams[k + nPVecs][2*kp + 1]);
				ptVecs[k] = curPtVec;
			}
		}
	}
	error |= discreteParams.size() != 3;
	if(!error)
	{
		VM::fileName.resize(discreteParams[0].size());
		for(std::size_t k = 0; k < VM::fileName.size(); ++k)
			VM::fileName[k] = discreteParams[0][k];
		if(discreteParams[1].size() == VM::vmTORSpec.size())
			VM::vmTORSpec = VolMeshTORSpecification(discreteParams[1]);
		else
			error = true;
		if(discreteParams[2].size() == VM::fixedBlockVec.size())
		{
			for(std::size_t k = 0; k < realParams[3].size(); ++k)
				VM::fixedBlockVec[k].first = discreteParams[2][k];
		}
		else
			error = true;
	}
	if(error)
	{
		std::cout << "Error in constructor for VolMesh2D, defining from a set of parameters" << std::endl;
		abort();
	}
	finishSetup(ptVecs);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::finishSetup(const std::vector<std::vector<Point_2_base> >& ptVecs)
{
	if(VM::vmTORSpec.torMeshType == mffPolygon)
	{
		// First set up segments
		boundaryVV.clear();
		if(!ptVecs.empty())
		{
			for(std::size_t k = 0; k < ptVecs.size(); ++k)
				GeometryTranslation::addPointsToSegments(ptVecs[k], boundaryVV);
		}
		// Set up bogus materials
		std::vector<GenericMaterial> matVec(1);
		// Mesh
		CDT_2 mesh = GeometryTranslation::mesh2D(boundaryVV, matVec, VM::meshParams);
		VM::upMesh = std::unique_ptr<TOMesh2D>(new TOMesh2D(std::move(mesh)));
		VM::fileName = "na";
	}
	else
	{
		if(VM::vmTORSpec.torMeshType == mffExodus)
			VM::upMesh = InputLoader::ExoTxtMeshLoader::loadExoTxtFileTri(VM::fileName);
		else if(VM::vmTORSpec.torMeshType == mffGMSH)
			VM::upMesh = InputLoader::GMSHTxtMeshLoader::loadGMSHTxtFile2D(VM::fileName);
		if(!VM::upMesh)
		{
			std::cout << "Error reading mesh: " << VM::fileName << ", aborting" << std::endl;
			abort();
		}
		findMeshEdges(); // Add mesh edges to boundaryVV
	}
	// Set up optimization parameters
	std::size_t numUnks = VM::vmTORSpec.torUnknownLocation == ulNode ? VM::upMesh->getNumNodes() : VM::upMesh->getNumElements();
	VM::pixelArray = std::vector<double>(numUnks, VM::threshold);
	// Set up filter
	VM::upFilt = std::unique_ptr<FilterBase>(new Filter2D<>(VM::upMesh.get(), VM::vmTORSpec.torUnknownLocation == ulElement));
	// Initialize fixed values
	VM::initFixedVals();
	// Set up filter derivative
	VM::computeDiffFilt();
}

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh2D<PenaltyFunc, ProjectionFunc>::~VolMesh2D()
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TopOptRep> VolMesh2D<PenaltyFunc, ProjectionFunc>::clone() const
{
	return std::unique_ptr<TopOptRep>(new VolMesh2D<PenaltyFunc, ProjectionFunc>(*this));
}

// Decode
template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::get2DSegments(std::vector<Mesh_Segment_2>& segVec) const
{
	// Compute element centroids
	std::vector<Point_2_base> ptVec = VM::getElemCentroids2D();
	// Add zeros out of bounds
	std::vector<std::vector<Mesh_Segment_2>> boundarySegs;
	findMeshEdges(boundarySegs);
	std::vector<double> zVals = VM::getProjElemDensities();
	VM::setFixedVals(zVals);
	for(std::size_t k1 = 0; k1 < boundarySegs.size(); ++k1)
	{
		const std::vector<Mesh_Segment_2>& curPoly = boundarySegs[k1];
		for(std::size_t k2 = 0; k2 < curPoly.size(); ++k2)
		{
			const Mesh_Segment_2& curSeg = curPoly[k2];
			Mesh_Vector_2 vSeg = curSeg.to_vector(), nSeg = vSeg.perpendicular(CGAL::NEGATIVE);
			double slen = sqrt(curSeg.squared_length()), nlen = sqrt(nSeg.squared_length());
			// Center of edge
			Point_2_base newPt = curSeg.source() + 0.5*vSeg + 0.5*(slen/nlen)*nSeg;
			zVals.push_back(0.);
			ptVec.push_back(newPt);
			// Point at the average of the normals of this edge and next
			std::size_t k2p1 = (k2 + 1) % curPoly.size();
			const Mesh_Segment_2& nextSeg = curPoly[k2p1];
			Mesh_Vector_2 vSeg2 = nextSeg.to_vector(), nSeg2 = vSeg2.perpendicular(CGAL::NEGATIVE);
			double slen2 = sqrt(nextSeg.squared_length()), nlen2 = sqrt(nSeg2.squared_length());
			Mesh_Vector_2 meanN = 0.5*(nSeg/nlen + nSeg2/nlen2);
			double nlenMean = sqrt(meanN.squared_length()), slenMean = 0.5*(slen + slen2);
			newPt = nextSeg.source() + 0.5*(slenMean/nlenMean)*meanN;
			zVals.push_back(0.);
			ptVec.push_back(newPt);
		}
	}
	segVec = PostProcess::meshSegsIsoSurf2d(zVals, ptVec, VM::threshold);
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TOMesh> VolMesh2D<PenaltyFunc, ProjectionFunc>::get2DMesh() const
{
	// Mesh is already formed, just need to assign values to triangles
	std::vector<double> tmpVec = VM::getProjElemDensities();
	VM::setFixedVals(tmpVec);
	// Penalize
	PenaltyFunc pfunc(VM::penalParams);
	std::transform(tmpVec.begin(), tmpVec.end(), tmpVec.begin(), pfunc);
	// Set new mesh
	std::unique_ptr<TOMesh> upOutMesh(VM::upMesh->clone());
	VM::setMeshOptVals(tmpVec, upOutMesh.get());
	return upOutMesh;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TOMesh> VolMesh2D<PenaltyFunc, ProjectionFunc>::get2DMesh(const GeometryTranslation::MesherData& meshParams) const
{
	return get2DMesh();
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TOMesh> VolMesh2D<PenaltyFunc, ProjectionFunc>::getOutputMesh() const
{
	std::vector<double> tmpVec = VM::getProjElemDensities();
	VM::setFixedVals(tmpVec);
	std::unique_ptr<TOMesh> upOutMesh(VM::upMesh->clone());
	VM::setMeshOptVals(tmpVec, upOutMesh.get());
	return upOutMesh;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::getBoundary(std::vector<Mesh_Segment_2>& boundaryVec) const
{
	boundaryVec.clear();
	GeometryTranslation::flattenMeshSegments(boundaryVV, boundaryVec);
}

// Modify
template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::refine()
{
	// Not implemented yet
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::prune()
{
	// TODO: Not implemented yet
	// No clear way to prune anyway
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::getDataSize(std::vector<std::size_t>& sizes) const
{
	sizes = {VM::pixelArray.size(), 1};
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
	std::vector<std::vector<double> >& realParams) const
{
	// Integer vals
	std::vector<int> fnameVec(VM::fileName.size());
	for(std::size_t k = 0; k < fnameVec.size(); ++k)
		fnameVec[k] = VM::fileName[k];
	discreteParams.resize(3);
  discreteParams[0] = fnameVec;
  discreteParams[1] = VM::vmTORSpec.toVec();
	discreteParams[2].resize(VM::fixedBlockVec.size());
	for(std::size_t k = 0; k < VM::fixedBlockVec.size(); ++k)
		discreteParams[2][k] = VM::fixedBlockVec[k].first;
	// Real vals
	std::size_t nPVecs = 4;
	realParams.resize(nPVecs + boundaryVV.size());
	realParams[0] = {VM::threshold, VM::filtRad, VM::meshParams.triMeshEdgeAngle, VM::meshParams.triMeshEdgeSize};
	realParams[1] = VM::penalParams;
	realParams[2] = VM::projParams;
	// Insert fixed block info
  realParams[3].resize(VM::fixedBlockVec.size());
  for(std::size_t k = 0; k < VM::fixedBlockVec.size(); ++k)
    realParams[3][k] = VM::fixedBlockVec[k].second;
	// Insert boundary points
	for(std::size_t k = 0; k < boundaryVV.size(); ++k)
	{
		std::vector<double> ptVec(2*boundaryVV[k].size());
		for(std::size_t kp = 0; kp < boundaryVV[k].size(); ++kp)
		{
			ptVec[2*kp] = boundaryVV[k][kp].source().x();
			ptVec[2*kp + 1] = boundaryVV[k][kp].source().y();
		}
		realParams[k + nPVecs] = ptVec;
	}
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::findMeshEdges()
{
	findMeshEdges(boundaryVV);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::findMeshEdges(std::vector<std::vector<Mesh_Segment_2>>& boundarySegs) const
{
	// Add mesh edges to boundaryVV
	// Search upMesh for boundary segments
	// Add each new edge to vector
	// If duplicates are found, this is an interior edge so delete it
	typedef std::pair<std::size_t, std::size_t> Edge_key;
	typedef std::unordered_set<Edge_key, Edge_hash, Edge_equal> Edge_set;
	Edge_set edgeConSet;
	for(std::size_t k = 0; k < VM::upMesh->getNumElements(); ++k)
	{
		std::vector<std::size_t> connVec = VM::upMesh->getElementConnectivity(k);
		for(std::size_t kn = 0; kn < connVec.size(); ++kn)
		{
			std::size_t knp1 = (kn + 1) % connVec.size();
			Edge_key curEdge(connVec[kn], connVec[knp1]);
			std::pair<Edge_set::iterator, bool> res = edgeConSet.insert(curEdge);
			if(!res.second)
			{
				// Edge has already been seen, remove it since it's an internal edge
				edgeConSet.erase(res.first);
			}
		}
	}
	// Get all points for the boundary edges
	std::vector<Mesh_Segment_2> edgeVec;
	for(Edge_set::const_iterator kit = edgeConSet.begin(); kit != edgeConSet.end(); ++kit)
	{
		Edge_key curEdge = *kit;
		edgeVec.push_back(Mesh_Segment_2(VM::upMesh->getNode2D(curEdge.first), VM::upMesh->getNode2D(curEdge.second)));
	}
	// Finally, merge boundary edges into boundarySegs
	GeometryTranslation::orderMeshSegments(edgeVec, boundarySegs);
	// Find outer boundary
	if(!boundarySegs.empty())
	{
		// First find bounding polygon with max area and save orientation of each
		std::size_t kmax = 0;
		double maxArea = GeometryTranslation::computeSignedArea(boundarySegs[0]);
		std::vector<bool> orientation(boundarySegs.size());
		orientation[0] = maxArea >= 0.;
		for(std::size_t k = 1; k < boundarySegs.size(); ++k)
		{
			double curArea = GeometryTranslation::computeSignedArea(boundarySegs[k]);
			orientation[k] = curArea >= 0.;
			if(fabs(curArea) > fabs(maxArea))
			{
				maxArea = curArea;
				kmax = k;
			}
		}
		// Correct orders
		for(std::size_t k = 0; k < boundarySegs.size(); ++k)
    {
			if((k == kmax && !orientation[k]) || (k != kmax && orientation[k]))
			{
				// Reverse vector and its segments
				// TODO: Put this in GeometryTranslation
				std::reverse(boundarySegs[k].begin(), boundarySegs[k].end());
				for(auto it = boundarySegs[k].begin(); it != boundarySegs[k].end(); ++it)
					*it = it->opposite();
			}
		}
	}
}

template class VolMesh2D<HelperNS::powPenal, HelperNS::defaultProjFunc>;
template class VolMesh2D<HelperNS::powPenalMin, HelperNS::defaultProjFunc>;
template class VolMesh2D<HelperNS::powPenal, HelperNS::regularizedHeaviside>;
template class VolMesh2D<HelperNS::powPenalMin, HelperNS::regularizedHeaviside>;
template class VolMesh2D<HelperNS::powPenal, HelperNS::thresholdHeaviside>;
template class VolMesh2D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>;
}
