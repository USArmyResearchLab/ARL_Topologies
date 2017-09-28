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

#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <string>
#include <vector>
#include "UTIL/topologiesdefs.h"
#include "REP/cgal_types.h"
#include "geometricentity.h"

namespace Topologies{
struct TOMesh;
struct TOMesh3D;
struct TOMesh2D;
}

//! Class that defines a boundary condition as a set of nodes from a mesh
class BoundaryCondition
{
public:
	//! Constructor defining a BoundaryCondition from an Exodus mesh file
	BoundaryCondition(BCType inBC, bool inFX, bool inFY, bool inFZ, unsigned nodeSetID, Topologies::MeshFileFormat inMFF,
		const std::string& meshFileName, unsigned dim);
	//! Constructor defining a 2D BoundaryCondition from a GeometricEntity file
	BoundaryCondition(BCType inBC, bool inFX, bool inFY, std::unique_ptr<GeometricEntity> inGE);
	//! Constructor defining a 3D BoundaryCondition from a GeometricEntity file
	BoundaryCondition(BCType inBC, bool inFX, bool inFY, bool inFZ, std::unique_ptr<GeometricEntity> inGE);
	BoundaryCondition(const BoundaryCondition& inBC);
	BoundaryCondition(BoundaryCondition&& rhs);
	BoundaryCondition& operator=(BoundaryCondition rhs);
	void swap(BoundaryCondition& rhs);
	~BoundaryCondition();
	//! Determines the set of nodes of inMesh corresponding to the boundary condition defined by this object
	std::vector<std::size_t> applyBC(const Topologies::TOMesh * const inMesh) const;
	//! Returns whether or not any node of inMesh coincides with this boundary condition
	bool checkValidity(const Topologies::TOMesh* const inMesh) const;
	//! Returns the boundary condition type BCType
	BCType getBCT();
	//! Returns the fixed coordinates
	std::vector<bool> getFixedCoords() const;
private:
	BCType type;
	bool fixX, fixY, fixZ;
	std::unique_ptr<GeometricEntity> upGE;
	std::vector<std::size_t> nodeIDVec;
	static const double tol;
};

inline
BCType BoundaryCondition::getBCT()
{
	return type;
}

#endif

