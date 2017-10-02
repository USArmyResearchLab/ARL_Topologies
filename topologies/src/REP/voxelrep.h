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

#ifndef VOXELREP_H
#define VOXELREP_H

#include "volmesh3d.h"
#include "inputloaderrep.h"
#include <memory>

namespace Topologies{
class GenericMaterial;

//! A topology representation that discretizes a rectangular solid region into voxels
/*! This is the most basic topology representation in 3d.  A region is broken into rectangular voxels
 *  each of which has an optimization parameter.  Note that in order for this to be effective
 *  sensitivity filtering is necessary.  Meshes can use either tet or hex elements.
 *  This is a specialization of VolMesh3D (which implements arbitrary 3d meshes) and is also parameterized
 *  by a penalization function and a projection function.
 */
template <typename PenaltyFunc = HelperNS::powPenalMin, typename ProjectionFunc = HelperNS::defaultProjFunc>
class VoxelRep : public VolMesh3D<PenaltyFunc, ProjectionFunc>
{
public:
	VoxelRep(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, const std::vector<double>& inPenalParams,
			const std::vector<double>& inProjParams);
	VoxelRep(TORType inTORT, const std::vector<std::vector<int> >& discreteParams, 
			const std::vector<std::vector<double>>& realParams);
	virtual ~VoxelRep();
	VoxelRep(const VoxelRep& copy);
	VoxelRep(VoxelRep&& copy) : VolMesh3D<PenaltyFunc, ProjectionFunc>(copy.myTORT) {swap(copy);}
	VoxelRep& operator=(VoxelRep rhs){swap(rhs); return *this;}
	void swap(VoxelRep& arg2);
	virtual std::unique_ptr<TopOptRep> clone() const;
	//! @name Decode functions
	//@{
	virtual std::unique_ptr<TOMesh> get3DSurfaceMesh() const;
	// TODO: Need 3d version of getBoundary
	//@}
	//! @name Structural modification functions
	//@{
	virtual void refine();
	virtual void prune();
	//@}
	//! @name Data access
	//@{
	virtual void getDataSize(std::vector<std::size_t>& sizes) const;
	virtual void getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
			std::vector<std::vector<double> >& realParams) const;
	virtual std::string getName() const {return getClassName();}
	//! Returns TopOptRep implementation's class name: voxel
	static std::string getClassName() {return "voxel";}
	//@}

private:
	void finishSetup();
	void refineElementUnknowns();
	void refineNodalUnknowns();

	typedef VolMesh<PenaltyFunc, ProjectionFunc> VM;
	typedef VolMesh3D<PenaltyFunc, ProjectionFunc> VM3D;

	double width, length, height; // Region physical size
	unsigned nx, ny, nz;
	MeshElementType myMET;
};

template <typename PenaltyFunc, typename ProjectionFunc>
VoxelRep<PenaltyFunc, ProjectionFunc>::VoxelRep(const VoxelRep<PenaltyFunc, ProjectionFunc>& copy) :
	VolMesh3D<PenaltyFunc, ProjectionFunc>(copy),
	nx(copy.nx),
  ny(copy.ny),
  nz(copy.nz),
  width(copy.width),
  length(copy.length),
  height(copy.height),
  myMET(copy.myMET)
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VoxelRep<PenaltyFunc, ProjectionFunc>::swap(VoxelRep<PenaltyFunc, ProjectionFunc>& arg2)
{
	VM3D::swap(arg2);
	std::swap(nx, arg2.nx);
	std::swap(ny, arg2.ny);
	std::swap(nz, arg2.nz);
	std::swap(width, arg2.width);
	std::swap(length, arg2.length);
	std::swap(height, arg2.height);
	std::swap(myMET, arg2.myMET);
}
}
#endif


