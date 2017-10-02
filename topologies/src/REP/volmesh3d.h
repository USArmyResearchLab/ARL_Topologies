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

#ifndef VOLMESH3DREP_H
#define VOLMESH3DREP_H

#include "volmesh.h"
#include "inputloaderrep.h"
#include "helper.h"
#include <memory>

namespace Topologies{
class GenericMaterial;
class TOMesh;
class TOMesh3D;

//! A topology representation that is an arbitrary 3d mesh
/*! 
 *  This class implements an arbitrary 3d volume mesh.  Each element contains an optimization
 *  parameter that can be used to scale material properties.  Sensitivity filter may be necessary.
 *  In addition, the class uses two template parameters which are a penalization function and a
 *  projection function.  The penalization function should penalize intermediate values such as
 *  a SIMP or mod-SIMP (mod-SIMP is the default).  The projection function uses both a density
 *  filter and the specified projection to aid optimization.  Default is nothing, but it can be
 *  a Heaviside projection.
 */
template <typename PenaltyFunc = HelperNS::powPenalMin, typename ProjectionFunc = HelperNS::defaultProjFunc>
class VolMesh3D : public VolMesh<PenaltyFunc, ProjectionFunc>
{
public:
	explicit VolMesh3D(TORType inTORT) : VolMesh<PenaltyFunc, ProjectionFunc>(inTORT) {}
	VolMesh3D(TORType inTORT, const InputLoader::TORGenericMesh& inputParams, const std::vector<double>& inPenalParams,
			const std::vector<double>& inProjParams);
	VolMesh3D(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, const std::vector<double>& inPenalParams,
			const std::vector<double>& inProjParams);
	VolMesh3D(TORType inTORT, const std::vector<std::vector<int> >& discreteParams, 
			const std::vector<std::vector<double> >& realParams);
	virtual ~VolMesh3D();
	VolMesh3D(const VolMesh3D& copy);
	VolMesh3D(VolMesh3D&& copy) : VolMesh<PenaltyFunc, ProjectionFunc>(copy.myTORT) {swap(copy);}
	VolMesh3D& operator=(VolMesh3D rhs) {swap(rhs); return *this;}
	void swap(VolMesh3D& arg2);
	virtual std::unique_ptr<TopOptRep> clone() const;

	//! @name Decode functions
	//@{
	virtual void get2DSegments(std::vector<Mesh_Segment_2>& segVec) const {}
	virtual std::unique_ptr<TOMesh> get2DMesh() const {return nullptr;}
	virtual std::unique_ptr<TOMesh> get2DMesh(const GeometryTranslation::MesherData& meshParams) const {return nullptr;}
	virtual void getBoundary(std::vector<Mesh_Segment_2>& boundaryVec) const {}
	virtual std::unique_ptr<TOMesh> get3DSurfaceMesh() const;
	virtual std::unique_ptr<TOMesh> get3DVolumeMesh() const;
	virtual std::unique_ptr<TOMesh> getOutputMesh() const;
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
	virtual unsigned getDimension() const;
	virtual void getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
		std::vector<std::vector<double> >& realParams) const;
	virtual std::string getName() const {return getClassName();}
	//! Returns TopOptRep implementation's class name: volmesh3d
	static std::string getClassName() {return "volmesh3d";}
	//@}

private:
	void finishSetup();
	void setOptVals(const std::vector<double>& optValVec, TOMesh3D* pTOM) const;
	Tr_GT::Sphere_3 getBoundingSphere(const TOMesh* const inMesh) const;
	std::vector<double> getElemAvgNodalDensities(const std::vector<double>& elemValVec) const;
	std::vector<Point_3_base> getElemCentroids() const;
	std::pair<double,double> getElementVolumeRange() const;
	void setFixedValsNodal(std::vector<double>& inVec) const;

	typedef VolMesh<PenaltyFunc, ProjectionFunc> VM;
protected:
	std::vector<double> getNodalDensities() const;

};

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh3D<PenaltyFunc, ProjectionFunc>::VolMesh3D(const VolMesh3D<PenaltyFunc, ProjectionFunc>& copy) :
	VolMesh<PenaltyFunc, ProjectionFunc>(copy)
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh3D<PenaltyFunc, ProjectionFunc>::swap(VolMesh3D<PenaltyFunc, ProjectionFunc>& arg2)
{
	VM::swap(arg2);
}

template <typename PenaltyFunc, typename ProjectionFunc>
unsigned VolMesh3D<PenaltyFunc, ProjectionFunc>::getDimension() const
{
	return 3;
}
}
#endif

