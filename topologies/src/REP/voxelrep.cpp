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

#include "voxelrep.h"
#include "geometrytranslation.h"
#include "helper.h"
#include "filter3d.h"
#include "tomesh.h"
#include "cartesianmesher.h"
#include "postprocess.h"
#include "cgal_types.h"
#include <algorithm>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

namespace Topologies{
template <typename PenaltyFunc, typename ProjectionFunc>
VoxelRep<PenaltyFunc, ProjectionFunc>::VoxelRep(TORType inTORT, const InputLoader::TORGenericVolume& inputParams,
		const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
	VolMesh3D<PenaltyFunc, ProjectionFunc>(inTORT, inputParams, inPenalParams, inProjParams),
	width(inputParams.getDimensions(0)),
	length(inputParams.getDimensions(1)),
	height(inputParams.getDimensions(2)),
	nx(inputParams.getDiscSizes(0)),
	ny(inputParams.getDiscSizes(1)),
	nz(inputParams.getDiscSizes(2)),
	myMET(inputParams.getMET())
{
	finishSetup();
}

template <typename PenaltyFunc, typename ProjectionFunc>
VoxelRep<PenaltyFunc, ProjectionFunc>::VoxelRep(TORType inTORT, const std::vector<std::vector<int> >& discreteParams, 
		const std::vector<std::vector<double> >& realParams) :
	VolMesh3D<PenaltyFunc, ProjectionFunc>(inTORT)
{
	bool error = discreteParams.empty();
	if(!error)
	{
		error = discreteParams[0].size() != 4;
		if(!error)
		{
			nx = discreteParams[0][0];
			ny = discreteParams[0][1];
			nz = discreteParams[0][2];
			myMET = (MeshElementType)discreteParams[0][3];
		}
		error |= discreteParams[1].size() != VM::vmTORSpec.size();
		if(!error)
			VM::vmTORSpec = VolMeshTORSpecification(discreteParams[1]);
	}
	error |= realParams.size() != 3;
	if(!error)
	{
		error = realParams[0].size() != 5;
		if(!error)
		{
			VM::threshold = realParams[0][0];
			VM::filtRad = realParams[0][1];
			width = realParams[0][2];
			length = realParams[0][3];
			height = realParams[0][4];
		}
		VM::penalParams = realParams[1];
		VM::projParams = realParams[2];
	}
	if(error)
	{
		std::cout << "Error in constructor VoxelRep, constructing from vectors of defining parameters." << std::endl;
		std::cout << "Constructor used for MPI code, check there for error!" << std::endl;
		abort();
	}
	finishSetup();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VoxelRep<PenaltyFunc, ProjectionFunc>::finishSetup()
{
	// Set up mesh
	std::size_t nelems = (std::size_t)nx*(std::size_t)ny*(std::size_t)nz;
	double initVal = 1./((double)nelems);
	std::vector<double> initArray = std::vector<double>(nelems, initVal);
	VM::upMesh = CartesianMesher::generateMesh(myMET, nx, ny, nz, width, length, height, initArray);
	// Set up unknown array
	if(VM::vmTORSpec.torUnknownLocation == ulNode)
		VM::pixelArray = std::vector<double>(VM::upMesh->getNumNodes(), initVal);
	else
		VM::pixelArray = initArray;
	// Set up filter
	VM::upFilt = std::unique_ptr<FilterBase>(new Filter3D<>(nx, ny, nz, width, length, height, 
		VM::vmTORSpec.torUnknownLocation == ulElement));
	// For tetrahedral meshes, there is not a 1-to-1 map between unknowns and elements
	// A filter is used to provide this functionality
	// Note that for highly skewed meshes (dx != dy != dz), this may be slightly incorrect as it will average close values
	if(myMET == metTet && VM::vmTORSpec.torUnknownLocation == ulElement && VM::filtRad == 0.)
	{
		VM::upDataFilt = std::move(VM::upFilt);
		VM::upFilt = std::unique_ptr<FilterBase>(new Filter3D<HelperNS::constantFunction>(nx, ny, nz, width, length, 
			height, VM::vmTORSpec.torUnknownLocation == ulElement));
		double dx = width/(double)nx, dy = length/(double)ny, dz = height/(double)nz;
//		VM::filtRad = (1. - 1e-6)*0.5*sqrt(dx*dx + dy*dy + dz*dz);
		VM::filtRad = 0.5*MAX(MAX(dx, dy), dz);
	}
	// Compute derivatives
	VM::computeDiffFilt();
}

template <typename PenaltyFunc, typename ProjectionFunc>
VoxelRep<PenaltyFunc, ProjectionFunc>::~VoxelRep()
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TopOptRep> VoxelRep<PenaltyFunc, ProjectionFunc>::clone() const
{
	return std::unique_ptr<TopOptRep>(new VoxelRep(*this));
}

// Modify
template <typename PenaltyFunc, typename ProjectionFunc>
void VoxelRep<PenaltyFunc, ProjectionFunc>::refine()
{
	// Refine grid, double linear sizes, octupling number of elements
	if(VM::vmTORSpec.torUnknownLocation == ulNode)
		refineNodalUnknowns();
	else
		refineElementUnknowns();
	nx *= 2;
	ny *= 2;
	nz *= 2;
	VM::upFilt = std::unique_ptr<FilterBase>(new Filter3D<>(nx, ny, nz, width, length, height, 
		VM::vmTORSpec.torUnknownLocation == ulElement));
	if(VM::upDataFilt)
	{
		VM::upDataFilt = std::move(VM::upFilt);
		VM::upFilt = std::unique_ptr<FilterBase>(new Filter3D<HelperNS::constantFunction>(nx, ny, nz, 
			width, length, height, VM::vmTORSpec.torUnknownLocation == ulElement));
		VM::filtRad *= 0.5;
	}
	VM::computeDiffFilt();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VoxelRep<PenaltyFunc, ProjectionFunc>::refineElementUnknowns()
{
	// New element optimization values are copied from old
	std::vector<double> newPixelArray(8*nx*ny*nz);
	std::size_t kxp = 0, kyp = 0, kzp = 0;
	for(std::size_t kz = 0; kz < 2*nz; ++kz)
	{
		for(std::size_t ky = 0; ky < 2*ny; ++ky)
		{
			for(std::size_t kx = 0; kx < 2*nx; ++kx)
			{
				newPixelArray[kz*4*nx*ny + ky*2*nx + kx] = VM::pixelArray[kzp*nx*ny + kyp*nx + kxp];
				kxp += kx % 2;
			}
			kxp = 0;
			kyp += ky % 2;
		}
		kyp = 0;
		kzp += kz % 2;
	}
	VM::pixelArray = newPixelArray;
	VM::upMesh = CartesianMesher::generateMesh(myMET, 2*nx, 2*ny, 2*nz, width, length, height, VM::pixelArray);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VoxelRep<PenaltyFunc, ProjectionFunc>::refineNodalUnknowns()
{
	// New element opt vals are interpolated using VoxelInterp
	std::vector<double> newPixelArray((2*nx + 1)*(2*ny + 1)*(2*nz + 1), 0.);
	VM::upMesh = CartesianMesher::generateMesh(myMET, 2*nx, 2*ny, 2*nz, width, length, height, 
										std::vector<double>(8*nx*ny*nz, 0.));
	PostProcess::VoxelInterp voxi(VM::pixelArray, nx, ny, nz, width, length, height, VM::threshold);
	for(std::size_t k = 0; k < VM::upMesh->getNumNodes(); ++k)
	{
		Point_3_base p = VM::upMesh->getNode3D(k);
		newPixelArray[k] = VM::threshold - voxi(p);
	}
	VM::pixelArray = newPixelArray;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VoxelRep<PenaltyFunc, ProjectionFunc>::prune()
{
	// TODO: Not implemented yet
	// Eventually will just reverse the refine process, not sure this is useful/necessary though
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::unique_ptr<TOMesh> VoxelRep<PenaltyFunc, ProjectionFunc>::get3DSurfaceMesh() const
{
  typedef CGAL::Implicit_surface_3<Tr_GT, PostProcess::VoxelInterp> TM_Surface_3;
  SM3_Tr tr;
  C2t3 c2t3(tr);
  PostProcess::VoxelInterp voxi(VM3D::getNodalDensities(), nx, ny, nz, width, length, height, VM::threshold);
	double rad2 = 0.25*(width*width + length*length + height*height);
	Tr_GT::Sphere_3 bounding_sphere(Point_3_base(width/2., length/2., height/2.), 1.5*rad2);
  TM_Surface_3 surface(voxi, bounding_sphere, 1e-5);
	double dx = width/(double)nx, dy = length/(double)ny, dz = height/(double)nz;
	double dxmin = MIN(MIN(dx, dy), dz);
  CGAL::Surface_mesh_default_criteria_3<SM3_Tr> criteria(30., 0.25*dxmin, 0.25*dxmin);
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
  return std::unique_ptr<TOMesh>(new TOMesh3DSurface(std::move(c2t3)));
}

// Data access
template <typename PenaltyFunc, typename ProjectionFunc>
void VoxelRep<PenaltyFunc, ProjectionFunc>::getDataSize(std::vector<std::size_t>& sizes) const
{
	sizes.resize(3);
	std::size_t addSize = VM::vmTORSpec.torUnknownLocation == ulNode ? 1 : 0;
	sizes[0] = nx + addSize;
	sizes[1] = ny + addSize;
	sizes[2] = nz + addSize;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VoxelRep<PenaltyFunc, ProjectionFunc>::getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
                                     std::vector<std::vector<double> >& realParams) const
{
	discreteParams.resize(2);
	discreteParams[0] = {(int)nx, (int)ny, (int)nz, (int)myMET};
	discreteParams[1] = VM::vmTORSpec.toVec();
	realParams.resize(3);
	realParams[0] = {VM::threshold, VM::filtRad, width, length, height};
	realParams[1] = VM::penalParams;
	realParams[2] = VM::projParams;
}

template class VoxelRep<HelperNS::powPenal, HelperNS::defaultProjFunc>;
template class VoxelRep<HelperNS::powPenalMin, HelperNS::defaultProjFunc>;
template class VoxelRep<HelperNS::powPenal, HelperNS::regularizedHeaviside>;
template class VoxelRep<HelperNS::powPenalMin, HelperNS::regularizedHeaviside>;
template class VoxelRep<HelperNS::powPenal, HelperNS::thresholdHeaviside>;
template class VoxelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>;
}
