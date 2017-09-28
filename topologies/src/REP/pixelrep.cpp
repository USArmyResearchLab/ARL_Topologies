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

#include "pixelrep.h"
#include "filter2d.h"
#include "tomesh.h"
#include "cartesianmesher.h"

namespace Topologies{
template <typename PenaltyFunc, typename ProjectionFunc>
PixelRep<PenaltyFunc, ProjectionFunc>::PixelRep(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, 
	const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
	VolMesh2D<PenaltyFunc, ProjectionFunc>(inTORT, inputParams, inPenalParams, inProjParams),
	width(inputParams.getDimensions(0)),
	height(inputParams.getDimensions(1)),
	nx(inputParams.getDiscSizes(0)),
	ny(inputParams.getDiscSizes(1)),
	myMET(inputParams.getMET())
{
	finishSetup();
}

template <typename PenaltyFunc, typename ProjectionFunc>
PixelRep<PenaltyFunc, ProjectionFunc>::PixelRep(TORType inTORT, const std::vector<std::vector<int> >& discreteParams, const std::vector<std::vector<double> >& realParams) :
	VolMesh2D<PenaltyFunc, ProjectionFunc>(inTORT)
{
	bool error = discreteParams.size() != 2;
	if(!error)
	{
		error = discreteParams[0].size() != 3;
		if(!error)
		{
			nx = discreteParams[0][0];
			ny = discreteParams[0][1];
			myMET = (MeshElementType)discreteParams[0][2];
		}
		error |= discreteParams[1].size() != VM::vmTORSpec.size();
		if(!error)
			VM::vmTORSpec = VolMeshTORSpecification(discreteParams[1]);
	}
	error |= realParams.size() != 3;
	if(!error)
	{
		error = realParams[0].size() != 4;
		if(!error)
		{
			VM::threshold = realParams[0][0];
			VM::filtRad = realParams[0][1];
			width = realParams[0][2];
			height = realParams[0][3];
		}
		VM::penalParams = realParams[1];
		VM::projParams = realParams[2];
	}
	if(error)
	{
		std::cout << "Error in constructor PixelRep, constructing from vectors of defining parameters." << std::endl;
		std::cout << "Constructor used for MPI code, check there for error!" << std::endl;
		abort();
	}
	finishSetup();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void PixelRep<PenaltyFunc, ProjectionFunc>::finishSetup()
{
	std::size_t nelems = (std::size_t)nx*(std::size_t)ny;
	double initVal = 1./((double)nelems);
	std::vector<double> initArray = std::vector<double>(nelems, initVal);
	VM::upMesh = CartesianMesher::generateMesh(myMET, nx, ny, width, height, initArray);
	if(VM::vmTORSpec.torUnknownLocation == ulElement)
		VM::pixelArray = initArray;
	else
		VM::pixelArray = std::vector<double>(VM::upMesh->getNumNodes(), initVal);
	VM::upFilt = std::unique_ptr<FilterBase>(new Filter2D<>(nx, ny, width, height, VM::vmTORSpec.torUnknownLocation == ulElement));
	// For triangular meshes, there is not a 1-to-1 map between unknowns and elements
  // A filter is used to provide this functionality
	// Note that for highly skewed meshes (dx != dy), this may be slightly incorrect as it will average close values
	if(myMET == metTri && VM::vmTORSpec.torUnknownLocation == ulElement && VM::filtRad == 0.)
	{
		VM::upDataFilt = std::move(VM::upFilt);	
		VM::upFilt = std::unique_ptr<FilterBase>(new Filter2D<HelperNS::constantFunction>(nx, ny, 
			width, height, VM::vmTORSpec.torUnknownLocation == ulElement));
		double dx = width/(double)nx, dy = height/(double)ny;
		VM::filtRad = 0.5*MAX(dx, dy);
	}
	// Set up boundaries
	VM2D::boundaryVV.resize(1);
	VM2D::boundaryVV[0] = {Mesh_Segment_2(Point_2_base(0., 0.), Point_2_base(width, 0.)),
									 Mesh_Segment_2(Point_2_base(width, 0.), Point_2_base(width, height)),
									 Mesh_Segment_2(Point_2_base(width, height), Point_2_base(0., height)),
									 Mesh_Segment_2(Point_2_base(0., height), Point_2_base(0., 0.))};
	// Compute Jacobian
	VM::computeDiffFilt();
}

template <typename PenaltyFunc, typename ProjectionFunc>
PixelRep<PenaltyFunc, ProjectionFunc>::~PixelRep()
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TopOptRep> PixelRep<PenaltyFunc, ProjectionFunc>::clone() const
{
	return std::unique_ptr<TopOptRep>(new PixelRep(*this));
}

// Modify
template <typename PenaltyFunc, typename ProjectionFunc>
void PixelRep<PenaltyFunc, ProjectionFunc>::refine()
{
	if(VM::vmTORSpec.torUnknownLocation == ulNode)
		refineNodalUnknowns();
	else
		refineElementUnknowns();
	nx *= 2;
	ny *= 2;
	VM::upMesh = CartesianMesher::generateMesh(myMET, nx, ny, width, height, std::vector<double>(nx*ny, 1.));
	VM::upFilt = std::unique_ptr<FilterBase>(new Filter2D<>(nx, ny, width, height, VM::vmTORSpec.torUnknownLocation == ulElement));
	if(VM::upDataFilt)
	{
		VM::upDataFilt = std::move(VM::upFilt);
		VM::upFilt = std::unique_ptr<FilterBase>(new Filter2D<HelperNS::constantFunction>(nx, ny,	width, height, 
			VM::vmTORSpec.torUnknownLocation == ulElement));
		VM::filtRad *= 0.5;
	}
	VM::computeDiffFilt();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void PixelRep<PenaltyFunc, ProjectionFunc>::refineElementUnknowns()
{
	// Refine grid, double linear sizes, quadruple number of elements
	// Densities are copied into new elements
	std::vector<double> newPixelArray(4*nx*ny);
	std::size_t kxp = 0, kyp = 0;
	for(std::size_t ky = 0; ky < 2*ny; ++ky)
	{
		for(std::size_t kx = 0; kx < 2*nx; ++kx)
		{
			newPixelArray[ky*2*nx + kx] = VM::pixelArray[kyp*nx + kxp];
			kxp += kx % 2;
		}
		kxp = 0;
		kyp += ky % 2;
	}
	VM::pixelArray = newPixelArray;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void PixelRep<PenaltyFunc, ProjectionFunc>::refineNodalUnknowns()
{
	// Refine grid, nodal unknowns must be averaged onto new nodes
	std::vector<double> newPixelArray((2*nx + 1)*(2*ny + 1));
	// Split up into 3 parts: new nodes along vertical lines,
	//  new nodes along horizontal lines, and new nodes within old elements
	// First, horizontal lines
	std::vector<double>& locPA = VM::pixelArray;
	for(std::size_t ky = 0; ky < ny + 1; ++ky)
		for(std::size_t kx = 0; kx < nx; ++kx)
			newPixelArray[2*ky*(2*nx+1) + 2*kx + 1] = 0.5*(locPA[(nx+1)*ky + kx] + locPA[(nx+1)*ky + kx + 1]);
	// Next, vertical lines
	for(std::size_t ky = 0; ky < ny; ++ky)
		for(std::size_t kx = 0; kx < nx + 1; ++kx)
			newPixelArray[(2*ky + 1)*(2*nx+1) + 2*kx] = 0.5*(locPA[(nx+1)*ky + kx] + locPA[(nx+1)*(ky + 1) + kx]);
	// Centers of old elements
	for(std::size_t ky = 0; ky < ny; ++ky)
    for(std::size_t kx = 0; kx < nx; ++kx)
      newPixelArray[(2*ky + 1)*(2*nx+1) + 2*kx + 1] = 0.25*(locPA[(nx+1)*ky + kx] + locPA[(nx+1)*ky + kx + 1] + 
																														locPA[(nx+1)*(ky + 1) + kx] + locPA[(nx+1)*(ky + 1) + kx + 1]);
	// Finally, copy old values to new nodes
	for(std::size_t ky = 0; ky < ny+1; ++ky)
		for(std::size_t kx = 0; kx < nx+1; ++kx)
			newPixelArray[(2*ky)*(2*nx+1) + 2*kx] = locPA[(nx+1)*ky + kx];
	locPA = newPixelArray;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void PixelRep<PenaltyFunc, ProjectionFunc>::prune()
{
	// TODO: Not implemented yet
	// Eventually will just reverse the refine process, not sure this is useful/necessary though
}

// Data access

template <typename PenaltyFunc, typename ProjectionFunc>
void PixelRep<PenaltyFunc, ProjectionFunc>::getDataSize(std::vector<std::size_t>& sizes) const
{
	sizes.resize(2);
	std::size_t addSize = VM::vmTORSpec.torUnknownLocation == ulNode ? 1 : 0;
	sizes[0] = nx + addSize;
	sizes[1] = ny + addSize;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void PixelRep<PenaltyFunc, ProjectionFunc>::getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
                                     std::vector<std::vector<double> >& realParams) const
{
	// Get defining parameters
	discreteParams.resize(2);
	discreteParams[0] = {(int)nx, (int)ny, (int)myMET};
	discreteParams[1] = VM::vmTORSpec.toVec();
	realParams.resize(3);
	realParams[0] = {VM::threshold, VM::filtRad, width, height};
	realParams[1] = VM::penalParams;
	realParams[2] = VM::projParams;
}

template class PixelRep<HelperNS::powPenal, HelperNS::defaultProjFunc>;
template class PixelRep<HelperNS::powPenalMin, HelperNS::defaultProjFunc>;
template class PixelRep<HelperNS::powPenal, HelperNS::regularizedHeaviside>;
template class PixelRep<HelperNS::powPenalMin, HelperNS::regularizedHeaviside>;
template class PixelRep<HelperNS::powPenal, HelperNS::thresholdHeaviside>;
template class PixelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>;
}
