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

#include "volmesh3d.h"
#include "geometrytranslation.h"
#include "helper.h"
#include "filter3d.h"
#include "tomesh.h"
#include "postprocess.h"
#include "outputwriter.h"
#include "exotxtmeshloader.h"
#include "gmshtxtloader.h"
#include "stltxtloader.h"
#include "tomeshprocessing.h"
//#include <random>
#include <algorithm>
#include <functional>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

namespace Topologies{
template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh3D<PenaltyFunc, ProjectionFunc>::VolMesh3D(TORType inTORT, const InputLoader::TORGenericVolume& inputParams,
    const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
  VolMesh<PenaltyFunc, ProjectionFunc>(inTORT, inputParams, inPenalParams, inProjParams)
{
	// Mesh initialization is handled by VoxelRep
}

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh3D<PenaltyFunc, ProjectionFunc>::VolMesh3D(TORType inTORT, const InputLoader::TORGenericMesh& inputParams,
		const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
	VolMesh<PenaltyFunc, ProjectionFunc>(inTORT, inputParams, inPenalParams, inProjParams)
{
	finishSetup();
}

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh3D<PenaltyFunc, ProjectionFunc>::VolMesh3D(TORType inTORT, const std::vector<std::vector<int> >& discreteParams, 
		const std::vector<std::vector<double> >& realParams) :
	VolMesh<PenaltyFunc, ProjectionFunc>(inTORT)
{
	bool error = discreteParams.size() != 3;
	if(!error)
	{
		VM::fileName.resize(discreteParams[0].size());
		for(std::size_t k = 0; k < VM::fileName.size(); ++k)
			VM::fileName[k] = discreteParams[0][k];
		if(discreteParams[1].size() == VM::vmTORSpec.size())
			VM::vmTORSpec = VolMeshTORSpecification(discreteParams[1]);
		else
			error = true;
		VM::fixedBlockVec.resize(discreteParams[2].size());
		for(std::size_t k = 0; k < VM::fixedBlockVec.size(); ++k)
			VM::fixedBlockVec[k].first = discreteParams[2][k];
	}
	error |= realParams.size() != 4;
	if(!error)
	{
		error = realParams[0].size() != 8;
		if(!error)
		{
			VM::threshold = realParams[0][0];
			VM::filtRad = realParams[0][1];
			VM::meshParams.tetMeshCellSize = realParams[0][2];
			VM::meshParams.tetMeshEdgeSize = realParams[0][3];
			VM::meshParams.tetMeshFacetAngle = realParams[0][4];
			VM::meshParams.tetMeshFacetSize = realParams[0][5];
			VM::meshParams.tetMeshFacetDistance = realParams[0][6];
			VM::meshParams.tetMeshCellRadiusEdgeRatio = realParams[0][7];
		}
		VM::penalParams = realParams[1];
		VM::projParams = realParams[2];
		if(realParams[3].size() == VM::fixedBlockVec.size())
		{
			for(std::size_t k = 0; k < VM::fixedBlockVec.size(); ++k)
				VM::fixedBlockVec[k].second = realParams[3][k];
		}
		else
			error = true;
	}
	if(error)
	{
		std::cout << "Error in constructor VolMesh3D, constructing from vectors of defining parameters." << std::endl;
		std::cout << "Constructor used for MPI code, check there for error!" << std::endl;
		abort();
	}
	finishSetup();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh3D<PenaltyFunc, ProjectionFunc>::finishSetup()
{
	// Set up mesh
	if(VM::vmTORSpec.torMeshType == mffExodus)
		VM::upMesh = InputLoader::ExoTxtMeshLoader::loadExoTxtFileTet(VM::fileName);
	else if(VM::vmTORSpec.torMeshType == mffGMSH)
		VM::upMesh = InputLoader::GMSHTxtMeshLoader::loadGMSHTxtFile3D(VM::fileName);
	else if(VM::vmTORSpec.torMeshType == mffSTL)
	{
		C3t3 cgalMesh = GeometryTranslation::mesh3D(InputLoader::STLTxtLoader::loadFileCGAL(VM::fileName), VM::meshParams);
		VM::upMesh = std::unique_ptr<TOMesh3D>(new TOMesh3D(cgalMesh));
	}
	if(!VM::upMesh)
	{
		std::cout << "Error reading mesh: " << VM::fileName << ", aborting" << std::endl;
		abort();
	}
	// Set up optimization value array
	std::size_t numUnks = VM::vmTORSpec.torUnknownLocation == ulNode ? VM::upMesh->getNumNodes() : VM::upMesh->getNumElements();
	VM::pixelArray = std::vector<double>(numUnks, VM::threshold);
	// Set up filter
	VM::upFilt = std::unique_ptr<FilterBase>(new Filter3D<>(VM::upMesh.get(), VM::vmTORSpec.torUnknownLocation == ulElement));
	// Set up fixed value blocks
	VM::initFixedVals();
	// Set up filter derivative, filter is linear so only need to do it once
	VM::computeDiffFilt();
}

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh3D<PenaltyFunc, ProjectionFunc>::~VolMesh3D()
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TopOptRep> VolMesh3D<PenaltyFunc, ProjectionFunc>::clone() const
{
	return std::unique_ptr<TopOptRep>(new VolMesh3D<PenaltyFunc, ProjectionFunc>(*this));
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh3D<PenaltyFunc, ProjectionFunc>::setFixedValsNodal(std::vector<double>& inVec) const
{
	if(VM::fixedBlockVec.empty())
		return;
	// Set up a vector of element-based fixed values
	std::vector<double> fixedVals(VM::pixelArray.size(), 0.);
	VM::setFixedVals(fixedVals);
	std::vector<double> fvmask = VM::getFixedValElemMask();
	// Interpolate those onto the nodes
	std::vector<double> nodalFixed = getElemAvgNodalDensities(fixedVals);
	std::vector<double> nodalMask = getElemAvgNodalDensities(fvmask);
	// Use mask as linear combination weight
	assert(inVec.size() == nodalFixed.size());
	assert(inVec.size() == nodalMask.size());
	for(std::size_t k = 0; k < inVec.size(); ++k)
		inVec[k] = nodalMask[k]*inVec[k] + (1. - nodalMask[k])*nodalFixed[k];
}

// Decode
template <typename PenaltyFunc, typename ProjectionFunc>
Tr_GT::Sphere_3 VolMesh3D<PenaltyFunc, ProjectionFunc>::getBoundingSphere(const TOMesh* const inMesh) const
{
	if(inMesh->getNumNodes() == 0)
	{
		std::cout << "Warning: Mesh has no nodes, in VolMesh3D<PenaltyFunc, ProjectionFunc>::getBoundingSphere()" << std::endl;
		return Tr_GT::Sphere_3(Tr_GT::Point_3(0., 0., 0.), 1.);
	}
	// Get radius
	Point_3_base p0 = inMesh->getNode3D(0);
	double minx = p0.x(), maxx = p0.x(), miny = p0.y(), maxy = p0.y(), minz = p0.z(), maxz = p0.z();
	for(std::size_t k = 1; k < inMesh->getNumNodes(); ++k)
	{
		p0 = inMesh->getNode3D(k);
		if(p0.x() < minx)
			minx = p0.x();
		else if(p0.x() > maxx)
			maxx = p0.x();
		if(p0.y() < miny)
			miny = p0.y();
		else if(p0.y() > maxy)
			maxy = p0.y();
		if(p0.z() < minz)
			minz = p0.z();
		else if(p0.z() > maxz)
			maxz = p0.z();
	}
	double rad2 = ((maxx - minx)*(maxx - minx) + (maxy - miny)*(maxy - miny) + (maxz - minz)*(maxz - minz));
	// Get center, must be inside the bounding surface
	std::size_t inElem = 0;
	bool found = false;
	for(std::size_t k = 0; k < inMesh->getNumElements() && !found; ++k)
	{
		if(inMesh->getOptVal(k) > VM::threshold)
		{
			inElem = k;
			found = true;
		}
	}
	Tr_GT::Point_3 center = TOMeshProcessing::getElementCentroid3D(inElem, inMesh);
	return Tr_GT::Sphere_3(center, 1.2*rad2); // 1.2 for safety factor
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh3D<PenaltyFunc, ProjectionFunc>::getNodalDensities() const
{
	if(VM::vmTORSpec.torUnknownLocation == ulNode)
	{
		// Nodal
		std::vector<double> res(VM::pixelArray.size());
		ProjectionFunc projf(VM::projParams);
		std::transform(VM::pixelArray.begin(), VM::pixelArray.end(), res.begin(), projf);
		setFixedValsNodal(res);
		return res;
	}
	// Elem-based
	std::vector<double> res = VM::getProjElemDensities();
	setFixedVals(res);
	std::vector<double> nodalVec = getElemAvgNodalDensities(res);	
	return nodalVec;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh3D<PenaltyFunc, ProjectionFunc>::getElemAvgNodalDensities(const std::vector<double>& elemValVec) const
{
	std::vector<double> nodalVec(VM::upMesh->getNumNodes(), 0.), volVec(VM::upMesh->getNumNodes(), 0.);
	for(std::size_t k = 0; k < VM::upMesh->getNumElements(); ++k)
	{
		// Take average of element-based densities
		std::vector<std::size_t> elemConn = VM::upMesh->getElementConnectivity(k);
		double vol = TOMeshProcessing::computeElementVolume(k, VM::upMesh.get());
		for(auto it = elemConn.begin(); it !=  elemConn.end(); ++it)
		{
			volVec[*it] += vol;
			nodalVec[*it] += vol*elemValVec[k];
		}
	}
	std::transform(nodalVec.begin(), nodalVec.end(), volVec.begin(), nodalVec.begin(), std::divides<double>());
	return nodalVec;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::pair<double,double> VolMesh3D<PenaltyFunc, ProjectionFunc>::getElementVolumeRange() const
{
	double minVol = 0., maxVol = 0.;
	if(VM::upMesh->getNumElements() == 0)
		return std::make_pair(minVol, maxVol);
	minVol = TOMeshProcessing::computeElementVolume(0, VM::upMesh.get());
	maxVol = minVol;
	for(std::size_t k = 1; k < VM::upMesh->getNumElements(); ++k)
	{
		double curVol = TOMeshProcessing::computeElementVolume(k, VM::upMesh.get());
		if(curVol < minVol)
			minVol = curVol;
		else if(curVol > maxVol)
			maxVol = curVol;
	}
	return std::make_pair(minVol, maxVol);
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TOMesh> VolMesh3D<PenaltyFunc, ProjectionFunc>::get3DSurfaceMesh() const
{
	// FEM-like interpolation on tetrahedra
	typedef CGAL::Implicit_surface_3<Tr_GT, PostProcess::TetMeshInterp> TM_Surface_3;
	SM3_Tr tr;
	C2t3 c2t3(tr);
	std::pair<double,double> volrange = getElementVolumeRange();
	double dxmin = pow(volrange.first*6.*sqrt(2.), 1./3.); // side len of a regular tet.
	double dxmax = pow(volrange.second*6.*sqrt(2.), 1./3.); // side len of a regular tet.
	std::unique_ptr<TOMesh> outMesh(getOutputMesh());
	PostProcess::TetMeshInterp toxi(getNodalDensities(), outMesh.get(), VM::threshold);
	Tr_GT::Sphere_3 bounding_sphere(getBoundingSphere(outMesh.get()));
	TM_Surface_3 surface(toxi, bounding_sphere, 1e-5);
	CGAL::Surface_mesh_default_criteria_3<SM3_Tr> criteria(30., 0.25*dxmin, 0.25*dxmin);
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
	std::unique_ptr<TOMesh> tmpTOMesh(new TOMesh3DSurface(std::move(c2t3)));
	return tmpTOMesh;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TOMesh> VolMesh3D<PenaltyFunc, ProjectionFunc>::get3DVolumeMesh() const
{
	// Mesh is already formed, just need to assign values to elements
	std::vector<double> tmpVec = VM::getProjElemDensities();
	VM::setFixedVals(tmpVec);
	// Penalty
	PenaltyFunc pfunc(VM::penalParams);
  std::transform(tmpVec.begin(), tmpVec.end(), tmpVec.begin(), pfunc);
	// Copy mesh
	std::unique_ptr<TOMesh> tmpTOMesh(VM::upMesh->clone());
	VM::setMeshOptVals(tmpVec, tmpTOMesh.get());
	return tmpTOMesh;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TOMesh> VolMesh3D<PenaltyFunc, ProjectionFunc>::getOutputMesh() const
{
	std::vector<double> tmpVec = VM::getProjElemDensities();
	VM::setFixedVals(tmpVec);
	std::unique_ptr<TOMesh> tmpTOMesh(VM::upMesh->clone());
	VM::setMeshOptVals(tmpVec, tmpTOMesh.get());
	return tmpTOMesh;
}

// Modify
template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh3D<PenaltyFunc, ProjectionFunc>::refine()
{
	// TODO: Not implemented yet
	// This will need to remesh the solid and interpolate the data onto the new nodes or elements
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh3D<PenaltyFunc, ProjectionFunc>::prune()
{
	// TODO: Not implemented yet
	// Eventually will just reverse the refine process, not sure this is useful/necessary though
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh3D<PenaltyFunc, ProjectionFunc>::getDataSize(std::vector<std::size_t>& sizes) const
{
	sizes = {VM::pixelArray.size(), 1, 1};
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh3D<PenaltyFunc, ProjectionFunc>::getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
	std::vector<std::vector<double> >& realParams) const
{
	// Get defining parameters
	std::vector<int> fnameVec(VM::fileName.size());
	for(std::size_t k = 0; k < fnameVec.size(); ++k)
		fnameVec[k] = VM::fileName[k];
	std::vector<int> ffVec = VM::vmTORSpec.toVec();
	std::vector<int> fvVec(VM::fixedBlockVec.size());
	for(std::size_t k = 0; k < VM::fixedBlockVec.size(); ++k)
    fvVec[k] = VM::fixedBlockVec[k].first;
	discreteParams = {fnameVec, ffVec, fvVec};
	// Real values
	realParams.resize(4);
	realParams[0] = {VM::threshold, VM::filtRad, VM::meshParams.tetMeshCellSize, VM::meshParams.tetMeshEdgeSize, 
									VM::meshParams.tetMeshFacetAngle, VM::meshParams.tetMeshFacetSize, VM::meshParams.tetMeshFacetDistance, 
									VM::meshParams.tetMeshCellRadiusEdgeRatio};
	realParams[1] = VM::penalParams;
	realParams[2] = VM::projParams;
	// Insert fixed block info
	realParams[3].resize(VM::fixedBlockVec.size());
	for(std::size_t k = 0; k < VM::fixedBlockVec.size(); ++k)
		realParams[3][k] = VM::fixedBlockVec[k].second;
}

template class VolMesh3D<HelperNS::powPenal, HelperNS::defaultProjFunc>;
template class VolMesh3D<HelperNS::powPenalMin, HelperNS::defaultProjFunc>;
template class VolMesh3D<HelperNS::powPenal, HelperNS::regularizedHeaviside>;
template class VolMesh3D<HelperNS::powPenalMin, HelperNS::regularizedHeaviside>;
template class VolMesh3D<HelperNS::powPenal, HelperNS::thresholdHeaviside>;
template class VolMesh3D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>;
}
