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

#ifndef VOLMESH2D_H
#define VOLMESH2D_H

#include "volmesh.h"
#include "inputloaderrep.h"
#include "helper.h"
#include <memory>
#include <functional>

namespace Topologies{
class GenericMaterial;
class TOMesh;
class TOMesh2D;

//! A topology representation that discretizes an arbitrary 2d region with tris or quads
/*! 
 *  This class implements an optimizable arbitrarily meshed 2d region.  Each element contains an
 *  optimization parameter that can be used to scale material properties.  Note that sensitivity
 *  filtering may be necessary.
 *  In addition, the class uses two template parameters which are a penalization function and a 
 *  projection function.  The penalization function should penalize intermediate values such as
 *  a SIMP or mod-SIMP (mod-SIMP is the default).  The projection function uses both a density
 *  filter and the specified projection to aid optimization.  Default is nothing, but it can be
 *  a Heaviside projection.
 */
template <typename PenaltyFunc = HelperNS::powPenalMin, typename ProjectionFunc = HelperNS::defaultProjFunc>
class VolMesh2D : public VolMesh<PenaltyFunc, ProjectionFunc>
{
public:
	explicit VolMesh2D(TORType inTORT) : VolMesh<PenaltyFunc, ProjectionFunc>(inTORT) {}
	VolMesh2D(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, const std::vector<double>& inPenalParams,
    const std::vector<double>& inProjParams);
	VolMesh2D(TORType inTORT, const InputLoader::TORGenericMesh& inputParams, const std::vector<double>& inPenalParams, 
		const std::vector<double>& inProjParams);
	VolMesh2D(TORType inTORT, const std::vector<std::vector<int>>& discreteParams, 
		const std::vector<std::vector<double> >& realParams);
	virtual ~VolMesh2D();
	VolMesh2D(const VolMesh2D& copy);
	VolMesh2D(VolMesh2D&& copy) : VolMesh<PenaltyFunc, ProjectionFunc>(copy.myTORT) {swap(copy);}
	VolMesh2D& operator=(VolMesh2D rhs) {swap(rhs); return *this;}
	void swap(VolMesh2D& arg2);
	virtual std::unique_ptr<TopOptRep> clone() const;

	//! @name Decode functions
	//@{
	virtual void get2DSegments(std::vector<Mesh_Segment_2>& segVec) const;
	virtual std::unique_ptr<TOMesh> get2DMesh() const;
	virtual std::unique_ptr<TOMesh> get2DMesh(const GeometryTranslation::MesherData& meshParams) const;
	virtual void getBoundary(std::vector<Mesh_Segment_2>& boundaryVec) const;
	virtual std::unique_ptr<TOMesh> getOutputMesh() const;
	//@}
	//! @name Structural modification functions
	//@{
	virtual void refine();
	virtual void prune();
	//@}
	//! @name Data access
	//@{
	virtual void getDataSize(std::vector<std::size_t>& sizes) const;
	virtual unsigned getDimension() const;
	virtual double getDataMagnitude() const {return 1.;}
	virtual void getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
		std::vector<std::vector<double> >& realParams) const;
	virtual std::string getName() const {return getClassName();}
	//! Returns TopOptRep implementation's class name: volmesh2d
	static std::string getClassName() {return "volmesh2d";}
	//@}

private:
	void finishSetup(const std::vector<std::vector<Point_2_base> >& ptVecs);
	void findMeshEdges();
	void findMeshEdges(std::vector<std::vector<Mesh_Segment_2>>& boundarySegs) const;

	typedef VolMesh<PenaltyFunc, ProjectionFunc> VM;
protected:
	std::vector<std::vector<Mesh_Segment_2>> boundaryVV;

private:
	// Hidden virtual functions
	virtual std::unique_ptr<TOMesh> get3DSurfaceMesh() const {return nullptr;}
	virtual std::unique_ptr<TOMesh> get3DVolumeMesh() const {return nullptr;}
};

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh2D<PenaltyFunc, ProjectionFunc>::VolMesh2D(const VolMesh2D<PenaltyFunc, ProjectionFunc>& copy) :
	VolMesh<PenaltyFunc, ProjectionFunc>(copy),
	boundaryVV(copy.boundaryVV)
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh2D<PenaltyFunc, ProjectionFunc>::swap(VolMesh2D<PenaltyFunc, ProjectionFunc>& arg2)
{
	VM::swap(arg2);
	boundaryVV.swap(arg2.boundaryVV);
}

template <typename PenaltyFunc, typename ProjectionFunc>
unsigned VolMesh2D<PenaltyFunc, ProjectionFunc>::getDimension() const
{
	return 2;
}

struct Edge_hash
{
	std::size_t operator()(std::pair<std::size_t, std::size_t> const& s) const
	{
		std::hash<std::size_t> uns_hash;
		return uns_hash(s.first) * uns_hash(s.second);
	}
};

struct Edge_equal
{
  bool operator()(std::pair<std::size_t, std::size_t> const& x, std::pair<std::size_t, std::size_t> const& y) const
  {
		return (x.first == y.first && x.second == y.second) || (x.first == y.second && x.second == y.first);
  }
};
}
#endif

