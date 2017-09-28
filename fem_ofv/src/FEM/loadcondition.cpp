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

#include "loadcondition.h"
#include "tomesh.h"
#include "exotxtmeshloader.h"
#include "gmshtxtloader.h"
#include <set>

// Hard-coded tolerance, maybe change to an input parameter
template <typename T>
const double LoadCondition<T>::tol = 1e-6;

using namespace Topologies;

template <typename T>
LoadCondition<T>::LoadCondition(BCType inLC, const std::vector<T>& inLoadVec, unsigned nodeSetID, MeshFileFormat inMFF, 
	const std::string& meshFileName, unsigned dim) :
	type(inLC),
  loadVec(inLoadVec)
{
	if(inMFF == mffExodus)
	{
		// Set up node IDs from the exodus mesh
		std::vector<std::pair<unsigned, std::vector<std::size_t>>> bcLoadVec;
		bcLoadVec = InputLoader::ExoTxtMeshLoader::loadNodeSets(meshFileName);
		if(bcLoadVec.empty())
		{
			std::cout << "Error: Didn't find any node sets in the supplied file: " << meshFileName << std::endl;
			abort();
		}
		for(std::size_t k = 0; k < bcLoadVec.size(); ++k)
		{
			if(nodeSetID == bcLoadVec[k].first)
				nodeIDVec = bcLoadVec[k].second;
		}
	}
	else if(inMFF == mffGMSH)
	{
		// Set up node IDs from the gmsh mesh
		std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> bcLoadVec;
		bcLoadVec = InputLoader::GMSHTxtMeshLoader::loadBoundaryPhysicalEntity(meshFileName, dim);
		if(bcLoadVec.empty())
		{
			std::cout << "Error: Didn't find any node sets in the supplied file: " << meshFileName << std::endl;
			abort();
		}
		for(std::size_t k = 0; k < bcLoadVec.size(); ++k)
		{
			if(nodeSetID == bcLoadVec[k].first)
			{
				// Add all node ids to set
				std::set<std::size_t> nodeSet;
				for(auto it1 = bcLoadVec[k].second.begin(); it1 != bcLoadVec[k].second.end(); ++it1)
					for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
						nodeSet.insert(*it2);
				nodeIDVec = std::vector<std::size_t>(nodeSet.begin(), nodeSet.end());
			}
		}
	}
	else
		std::cout << "Warning: LoadCondition received unsupported file type" << std::endl;
}

template <typename T>
LoadCondition<T>::LoadCondition(BCType inBC, const std::vector<T>& inLoadVec, std::unique_ptr<GeometricEntity> inGE):
	type(inBC),
	loadVec(inLoadVec),
	upGE(std::move(inGE))
{
}

template <typename T>
LoadCondition<T>::LoadCondition(const LoadCondition<T>& inLC):
	type(inLC.type),
	loadVec(inLC.loadVec),
	nodeIDVec(inLC.nodeIDVec)
{
	if(inLC.upGE)
		upGE = inLC.upGE->clone();
}

template <typename T>
LoadCondition<T>::LoadCondition(LoadCondition<T>&& rhs)
{
  swap(rhs);
}

template <typename T>
LoadCondition<T>& LoadCondition<T>::operator=(LoadCondition<T> rhs)
{
  swap(rhs);
  return *this;
}

template <typename T>
void LoadCondition<T>::swap(LoadCondition<T>& rhs)
{
  std::swap(loadVec, rhs.loadVec);
  std::swap(type, rhs.type);
  upGE.swap(rhs.upGE);
	nodeIDVec.swap(rhs.nodeIDVec);
}

template <typename T>
LoadCondition<T>::~LoadCondition()
{
}

template <typename T>
void LoadCondition<T>::applyLC(const TOMesh* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
                std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const
{
	if(inMesh->dimNum() == 2)
		applyLC(dynamic_cast<const TOMesh2D* const>(inMesh), elemIds, lcVecX, lcVecY);
	else if(inMesh->dimNum() == 3)
		applyLC(dynamic_cast<const TOMesh3D* const>(inMesh), elemIds, lcVecX, lcVecY, lcVecZ);
}

template <typename T>
void LoadCondition<T>::applyLC(const TOMesh2D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
                std::vector<T>& lcVecX, std::vector<T>& lcVecY) const
{
	if(upGE)
	{
		if(upGE->getGeomDim() == 0) // Add points
			applyLC0(inMesh, elemIds, lcVecX, lcVecY);
		else if(upGE->getGeomDim() == 1) // Add edges
			applyLC1(inMesh, elemIds, lcVecX, lcVecY);
	}
	// Add node sets
	lcVecX.insert(lcVecX.end(), nodeIDVec.size(), loadVec[0]);
	lcVecY.insert(lcVecY.end(), nodeIDVec.size(), loadVec[1]);
	for(std::size_t k = 0; k < nodeIDVec.size(); ++k)
	{
		std::vector<std::size_t> tmp(1, nodeIDVec[k]);
		elemIds.push_back(tmp);
	}
}

template <typename T>
void LoadCondition<T>::applyLC0(const TOMesh2D* const inMesh, std::vector< std::vector<std::size_t> >& elemIds,
                std::vector<T>& lcVecX, std::vector<T>& lcVecY) const
{
	for(std::size_t k = 0; k < inMesh->getNumNodes(); ++k)
	{
		Mesh_K::Point_2 tmpPt1 = inMesh->getNode2D(k);
		if(upGE->isPointCoincident(tmpPt1, tol))
		{
			std::vector<std::size_t> tmp;
			tmp.push_back(k);
			elemIds.push_back(tmp);
			lcVecX.push_back(loadVec[0]);
			lcVecY.push_back(loadVec[1]);
		}
	}
}

bool findEdge(const std::vector<std::vector<std::size_t>>& elemIds, const std::vector<std::size_t>& edge)
{
	if(edge.size() != 2)
		return false;
	for(std::size_t k = 0; k < elemIds.size(); ++k)
	{
		if(elemIds[k].size() != 2)
			return false;
		if((elemIds[k][0] == edge[0] && elemIds[k][1] == edge[1]) || (elemIds[k][0] == edge[1] && elemIds[k][1] == edge[0]))
			return true;
	}
	return false;
}

template <typename T>
void LoadCondition<T>::applyLC1(const TOMesh2D* const inMesh, std::vector<std::vector<std::size_t>>& elemIds, 
														std::vector<T>& lcVecX, std::vector<T>& lcVecY) const
{
	for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
	{
		std::vector<std::size_t> curElem = inMesh->getElementConnectivity(k);
		std::size_t numNodes = curElem.size();
		for(std::size_t ke = 0; ke < numNodes; ++ke)
		{
			std::size_t kep1 = (ke + 1)%numNodes;
			std::size_t fid1 = curElem[ke], fid2 = curElem[kep1];
			Mesh_K::Point_2 tmpPt1 = inMesh->getNode2D(fid1), tmpPt2 = inMesh->getNode2D(fid2);
			std::vector<std::size_t> tmp;
			if(upGE->isPointCoincident(tmpPt1, tol) && upGE->isPointCoincident(tmpPt2, tol))
			{
				tmp.push_back(fid1);
				tmp.push_back(fid2);
			}
			if(tmp.size() == 2)
			{
				if(!findEdge(elemIds, tmp))
					elemIds.push_back(tmp);
				lcVecX.push_back(loadVec[0]);
				lcVecY.push_back(loadVec[1]);
			}
		}
	}
}

template <typename T>
void LoadCondition<T>::applyLC(const TOMesh3D* const inMesh, std::vector<std::vector<std::size_t>>& elemIds, 
														std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const
{
	if(upGE)
	{
		if(upGE->getGeomDim() == 0) // Add points
			applyLC0(inMesh, elemIds, lcVecX, lcVecY, lcVecZ);
		else if(upGE->getGeomDim() == 1) // Add edges
			applyLC1(inMesh, elemIds, lcVecX, lcVecY, lcVecZ);
		else if(upGE->getGeomDim() == 2) // Add faces
			applyLC2(inMesh, elemIds, lcVecX, lcVecY, lcVecZ);
	}
	// Add node sets
	lcVecX.insert(lcVecX.end(), nodeIDVec.size(), loadVec[0]);
	lcVecY.insert(lcVecY.end(), nodeIDVec.size(), loadVec[1]);
	lcVecZ.insert(lcVecZ.end(), nodeIDVec.size(), loadVec[2]);
	for(std::size_t k = 0; k < nodeIDVec.size(); ++k)
	{
		std::vector<std::size_t> tmp(1, nodeIDVec[k]);
		elemIds.push_back(tmp);
	}
}

template <typename T>
void LoadCondition<T>::applyLC0(const TOMesh3D* const inMesh, std::vector<std::vector<std::size_t>>& elemIds,
                            std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const
{
	for(std::size_t k = 0; k < inMesh->getNumNodes(); ++k)
	{
		Mesh_K::Point_3 tmpP = inMesh->getNode3D(k);
		if(upGE->isPointCoincident(tmpP, tol))
		{
			std::vector<std::size_t> tmp;
			tmp.push_back(k);
			elemIds.push_back(tmp);
			lcVecX.push_back(loadVec[0]);
			lcVecY.push_back(loadVec[1]);
			lcVecZ.push_back(loadVec[2]);
		}
	}
}

template <typename T>
void LoadCondition<T>::applyLC1(const TOMesh3D* const inMesh, std::vector<std::vector<std::size_t>>& elemIds,
                            std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const
{
	for(std::size_t k = 0; k < inMesh->edgeVec.size(); ++k)
	{
		const std::vector<std::size_t>& edgeVec = inMesh->edgeVec[k];
		std::vector<std::size_t> tmp;
		for(unsigned k2 = 0; k2 < 2; k2++)
		{
			std::size_t nodeid = edgeVec[k2];
			Mesh_K::Point_3 tmpP = inMesh->getNode3D(nodeid);
			if(upGE->isPointCoincident(tmpP, tol))
				tmp.push_back(nodeid);
		}
		if(tmp.size() == 2)
		{
			elemIds.push_back(tmp);
			lcVecX.push_back(loadVec[0]);
			lcVecY.push_back(loadVec[1]);
			lcVecZ.push_back(loadVec[2]);
    }
  }
	if(lcVecX.empty())
		applyLC0(inMesh, elemIds, lcVecX, lcVecY, lcVecZ); // Found nothing, try nodes
}

template <typename T>
void LoadCondition<T>::applyLC2(const TOMesh3D* const inMesh, std::vector<std::vector<std::size_t>>& elemIds,
                            std::vector<T>& lcVecX, std::vector<T>& lcVecY, std::vector<T>& lcVecZ) const
{
	for(std::size_t k = 0; k < inMesh->faceVec.size(); ++k)
	{
		const std::vector<std::size_t>& faceVec = inMesh->faceVec[k];
		std::vector<std::size_t> tmp;
		for(unsigned k2 = 0; k2 < 3; k2++)
		{
			std::size_t nodeid = faceVec[k2];
			Mesh_K::Point_3 tmpP = inMesh->getNode3D(nodeid);
			if(upGE->isPointCoincident(tmpP, tol))
				tmp.push_back(nodeid);
		}
		if(tmp.size() == 3)
		{
			elemIds.push_back(tmp);
			lcVecX.push_back(loadVec[0]);
			lcVecY.push_back(loadVec[1]);
			lcVecZ.push_back(loadVec[2]);
		}
	}
	if(lcVecX.empty())
		applyLC1(inMesh, elemIds, lcVecX, lcVecY, lcVecZ); // Found nothing, try edges
}

template <typename T>
void LoadCondition<T>::applyLC(const TOMesh* const inMesh, std::vector<std::size_t>& nodeIds) const
{
	if(inMesh->dimNum() == 2)
		applyLC(dynamic_cast<const TOMesh2D* const>(inMesh), nodeIds);
	else if(inMesh->dimNum() == 3)
		applyLC(dynamic_cast<const TOMesh3D* const>(inMesh), nodeIds);
}

template <typename T>
void LoadCondition<T>::applyLC(const TOMesh2D* const inMesh, std::vector<std::size_t>& nodeIds) const
{
	if(upGE)
	{
		std::vector<std::vector<std::size_t>> elemIds;
		std::vector<T> lcVecX, lcVecY;
		applyLC0(inMesh, elemIds, lcVecX, lcVecY);
		nodeIds.resize(elemIds.size());
		for(std::size_t k = 0; k < elemIds.size(); ++k)
		{
			if(!elemIds[k].empty())
				nodeIds[k] = elemIds[k][0];
			else
				std::cout << "Warning! in applyLC, elemIds contained an empty element." << std::endl;
		}
	}
	// Add any nodeIDs
	nodeIds.insert(nodeIds.end(), nodeIDVec.begin(), nodeIDVec.end());
}

template <typename T>
void LoadCondition<T>::applyLC(const TOMesh3D* const inMesh, std::vector<std::size_t>& nodeIds) const
{
	if(upGE)
	{
		std::vector<std::vector<std::size_t>> elemIds;
		std::vector<T> lcVecX, lcVecY, lcVecZ;
		applyLC0(inMesh, elemIds, lcVecX, lcVecY, lcVecZ);
		nodeIds.resize(elemIds.size());
		for(std::size_t k = 0; k < elemIds.size(); ++k)
		{
			if(!elemIds[k].empty())
				nodeIds[k] = elemIds[k][0];
			else
				std::cout << "Warning! in applyLC, elemIds contained an empty element." << std::endl;
		}
	}
	// Add any nodeIDs
	nodeIds.insert(nodeIds.end(), nodeIDVec.begin(), nodeIDVec.end());
}

template <typename T>
bool LoadCondition<T>::checkValidity(const TOMesh* const inMesh) const
{
	if(!inMesh)
		return false;
	if(inMesh->dimNum() == 2)
		return checkValidity(dynamic_cast<const TOMesh2D* const>(inMesh));
	else if(inMesh->dimNum() == 3)
		return checkValidity(dynamic_cast<const TOMesh3D* const>(inMesh));
	return false;
}

template <typename T>
bool LoadCondition<T>::checkValidity(const TOMesh2D* const inMesh) const
{
	std::vector<std::vector<std::size_t>> elemIds;
	std::vector<T> lcVecX, lcVecY;
	applyLC(inMesh, elemIds, lcVecX, lcVecY);
	return !elemIds.empty();
}

template <typename T>
bool LoadCondition<T>::checkValidity(const TOMesh3D* const inMesh) const
{
	std::vector<std::vector<std::size_t>> elemIds;
	std::vector<T> lcVecX, lcVecY, lcVecZ;
	applyLC(inMesh, elemIds, lcVecX, lcVecY, lcVecZ);
	return !elemIds.empty();
}

template class LoadCondition<int>;
template class LoadCondition<double>;
