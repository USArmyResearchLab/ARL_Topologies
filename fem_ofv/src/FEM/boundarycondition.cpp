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

#include "boundarycondition.h"
#include "REP/tomesh.h"
#include "IO/exotxtmeshloader.h"
#include "IO/gmshtxtloader.h"

const double BoundaryCondition::tol = 1e-6;

using namespace Topologies;

BoundaryCondition::BoundaryCondition(BCType inBC, bool inFX, bool inFY, bool inFZ, unsigned nodeSetID, MeshFileFormat inMFF,
	const std::string& meshFileName, unsigned dim) :
	type(inBC),
	fixX(inFX),
	fixY(inFY),
	fixZ(inFZ)
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
		std::cout << "Warning: BoundaryCondition received unsupported file type" << std::endl;
}

BoundaryCondition::BoundaryCondition(BCType inBC, bool inFX, bool inFY, std::unique_ptr<GeometricEntity> inGE):
	type(inBC),
	fixX(inFX),
	fixY(inFY),
	fixZ(false),
	upGE(std::move(inGE))
{
}

BoundaryCondition::BoundaryCondition(BCType inBC, bool inFX, bool inFY, bool inFZ, std::unique_ptr<GeometricEntity> inGE):
	type(inBC),
	fixX(inFX),
	fixY(inFY),
	fixZ(inFZ),
	upGE(std::move(inGE))
{
}

BoundaryCondition::BoundaryCondition(const BoundaryCondition& inBC) :
	type(inBC.type),
	fixX(inBC.fixX),
	fixY(inBC.fixY),
	fixZ(inBC.fixZ),
	nodeIDVec(inBC.nodeIDVec)
{
	if(inBC.upGE)
		upGE = inBC.upGE->clone();
}

BoundaryCondition::BoundaryCondition(BoundaryCondition&& rhs)
{
	swap(rhs);
}

BoundaryCondition& BoundaryCondition::operator=(BoundaryCondition rhs)
{
	swap(rhs);
	return *this;
}

void BoundaryCondition::swap(BoundaryCondition& rhs)
{
	std::swap(fixX, rhs.fixX);
	std::swap(fixY, rhs.fixY);
	std::swap(fixZ, rhs.fixZ);
	std::swap(type, rhs.type);
	upGE.swap(rhs.upGE);
	nodeIDVec.swap(rhs.nodeIDVec);
}

BoundaryCondition::~BoundaryCondition()
{
}

std::vector<bool> BoundaryCondition::getFixedCoords() const
{
	std::vector<bool> outfix(3);
	outfix[0] = fixX;
	outfix[1] = fixY;
	outfix[2] = fixZ;
	return outfix;
}

std::vector<std::size_t> BoundaryCondition::applyBC(const TOMesh * const inMesh) const
{
	std::vector<std::size_t> outvec = nodeIDVec;
	if(upGE)
	{
		for(std::size_t kv = 0; kv < inMesh->getNumNodes(); ++kv)
		{
			if(inMesh->dimNum() == 2)
			{
				if(upGE->isPointCoincident(inMesh->getNode2D(kv), tol))
					outvec.push_back(kv);
			}
			else
			{
				if(upGE->isPointCoincident(inMesh->getNode3D(kv), tol))
					outvec.push_back(kv);
			}
		}
	}
	return outvec;
}

bool BoundaryCondition::checkValidity(const TOMesh* const inMesh) const
{
	if(!inMesh)
		return false;
	std::vector<std::size_t> tmp = applyBC(inMesh);
	return !tmp.empty();
}


