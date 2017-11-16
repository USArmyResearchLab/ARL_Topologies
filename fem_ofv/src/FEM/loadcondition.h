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

#ifndef LOADCONDITION_H
#define LOADCONDITION_H

#include <vector>
#include <map>
#include <string>
#include "topologiesdefs.h"
#include "REP/cgal_types.h"
#include "geometricentity.h"
#include "coordinatesystem.h"

namespace Topologies{
class TOMesh;
class TOMesh3D;
class TOMesh2D;
}

//! Class defining a loading condition as a set of nodes on a mesh
template<typename T>
class LoadCondition
{
public:
	//! Constructor defining a LoadCondition from an Exodus mesh file
	LoadCondition(BCType inLC, const std::vector<T>& inLoadVec, unsigned nodeSetID, Topologies::MeshFileFormat inMFF, 
		const std::string& meshFileName, unsigned dim, CoordinateSystem::Type inCT = CoordinateSystem::Type::cartesian);
	//! Constructor defining a LoadCondition from a GeometricEntity file
	LoadCondition(BCType inLC, const std::vector<T>& inLoadVec, std::unique_ptr<GeometricEntity> inGE, 
								CoordinateSystem::Type inCT = CoordinateSystem::Type::cartesian);
	LoadCondition(const LoadCondition& inLC);
  LoadCondition(LoadCondition&& rhs);
  LoadCondition& operator=(LoadCondition rhs);
  void swap(LoadCondition& rhs);
	~LoadCondition();
	//! Generates a set of elements and load vectors to apply
	/*! @param elemIds On return, contains a vector of elements, the size of each element is set by the dimensionality of the BC
	 *  @param lcVecX On return, contains the x-component of the load vector that corresponds th each element in elemIds
	 *  @param lcVecY On return, contains the y-component of the load vector that corresponds th each element in elemIds
	 *  @param lcVecZ On return, contains the z-component of the load vector that corresponds th each element in elemIds
	 */
	void applyLC(const Topologies::TOMesh* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
								std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const;
	//! Generates a set of nodes that the load is applied to
	void applyLC(const Topologies::TOMesh* const inMesh, std::vector<std::size_t>& nodeIds) const;
	//! Generates a set of elements and load vectors to apply
  /*! @param elemIds On return, contains a vector of elements, the size of each element is set by the dimensionality of the BC
	 *  @param lcVecX On return, contains the x-component of the load vector that corresponds th each element in elemIds
	 *  @param lcVecY On return, contains the y-component of the load vector that corresponds th each element in elemIds
	 */
	void applyLC(const Topologies::TOMesh2D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds, 
							std::vector<T>& lcVecX, std::vector<T>& lcVecY) const;
	//! Generates a set of nodes that the load is applied to
	void applyLC(const Topologies::TOMesh2D* const inMesh, std::vector<std::size_t>& nodeIds) const;
	//! Generates a set of elements and load vectors to apply
  /*! @param elemIds On return, contains a vector of elements, the size of each element is set by the dimensionality of the BC
	 *  @param lcVecX On return, contains the x-component of the load vector that corresponds th each element in elemIds
	 *  @param lcVecY On return, contains the y-component of the load vector that corresponds th each element in elemIds
	 *  @param lcVecZ On return, contains the z-component of the load vector that corresponds th each element in elemIds
	 */
	void applyLC(const Topologies::TOMesh3D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds, 
							std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const;
	//! Generates a set of nodes that the load is applied to
	void applyLC(const Topologies::TOMesh3D* const inMesh, std::vector<std::size_t>& nodeIds) const;
	//! Checks wether or not any element or node in inMesh matches this LoadCondition
	bool checkValidity(const Topologies::TOMesh* const inMesh) const;
	//! Checks wether or not any element or node in inMesh matches this LoadCondition
	bool checkValidity(const Topologies::TOMesh2D* const inMesh) const;
	//! Checks wether or not any element or node in inMesh matches this LoadCondition
	bool checkValidity(const Topologies::TOMesh3D* const inMesh) const;
	//! Returns the boundary condition type
	BCType getLCT() const;
	//! Returns the load vector
	std::vector<T> getLoadVec() const;
	//! Returns the coordinate system type of the load vector
	CoordinateSystem::Type getCoordinateSystemType() const;
private:
	void applyLC0(const Topologies::TOMesh2D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
								std::vector<T>& lcVecX, std::vector<T>& lcVecY) const;
	void applyLC1(const Topologies::TOMesh2D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
								std::vector<T>& lcVecX, std::vector<T>& lcVecY) const;
	void applyLC0(const Topologies::TOMesh3D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
                std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const;
	void applyLC1(const Topologies::TOMesh3D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
                std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const;
	void applyLC2(const Topologies::TOMesh3D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
                std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const;
private:
	BCType type;
	std::unique_ptr<GeometricEntity> upGE;
	std::vector<T> loadVec;
	std::vector<std::size_t> nodeIDVec;
	CoordinateSystem::Type ct;
	static const double tol;
};

template <typename T>
BCType LoadCondition<T>::getLCT() const
{
	return type;
}

template <typename T>
std::vector<T> LoadCondition<T>::getLoadVec() const
{
	return loadVec;
}

template <typename T>
CoordinateSystem::Type LoadCondition<T>::getCoordinateSystemType() const
{
	return ct;
}

#endif

