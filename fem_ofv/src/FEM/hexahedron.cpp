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

#include "hexahedron.h"
#include "linquad.h"
#include "lintetra.h"

Hexahedron::Hexahedron(CellType inCellType, const std::vector<Point3D*>& inNodeVec, const Topologies::GenericMaterial& inMat):
	Cell(inCellType, inMat)
{
	assert(inNodeVec.size() == 8);
	nodeVec = inNodeVec;
	setInitialNodeOrdering();
	patchVec = std::vector<Element<Point3D>*>(6, nullptr);
	bfVec.resize(nodeVec.size());
	for(std::size_t i = 0; i < bfVec.size(); i++)
		bfVec[i] = 0;
}

std::vector<std::size_t> Hexahedron::getNeighborNodes(std::size_t knode) const
{
	assert(knode < 8);
	// For a given node (rootNode), finds the closest neighbors (connectivity-wise) by 
	// finding the set of 3 nodes with maximal solid angle.  
	// Note this may not work for wonky shaped hexes, though you shouldn't be using them!
	// First, fill out a vector of node IDs to check
	std::vector<std::size_t> idVec;
	idVec.reserve(7);
	for(std::size_t k = 0; k < 8; ++k)
	{
		if(k != knode)
			idVec.push_back(k);
	}
	// Next, search all combinations (7 choose 3 = 35) of nodes for max solid angle triplet
	std::vector<bool> permVec(7,false);
	std::vector<std::size_t> curPts, maxSAPts;
	permVec[0] = permVec[1] = permVec[2] = true;
	double maxSolidAngle = 0.;
	const Point3D* chkPt = nodeVec[knode];
	curPts.reserve(3);
	do
	{
		for(std::size_t k = 0; k < 7; ++k)
		{
			if(permVec[k])
				curPts.push_back(idVec[k]);
		}
		assert(curPts.size() == 3);
		// Check for max solid angle
		double curSA = solidAngle(chkPt, *nodeVec[curPts[0]], *nodeVec[curPts[1]], *nodeVec[curPts[2]]);
		if(curSA > maxSolidAngle)
		{
			maxSolidAngle = curSA;
			maxSAPts = curPts;
		}
		curPts.clear();
	} while(prev_permutation(permVec.begin(), permVec.end()));
	return maxSAPts;
}

void Hexahedron::setInitialNodeOrdering()
{
	// Since a hexahedron is not a simplex, we need some pre-processing of the nodes to ensure correct ordering
	// Several functions below rely on this ordering, eg. createPatch
	// Order is something like the following, though it may be rotated:
	// 4*-------*7
	//  |\      |\
	//  |5*-------*6
	//  | |     | |
	// 0*-|-----*3|
	//   \|      \|
	//   1*-------*2
	// Essentially, one face should contain the first 4 nodes in the same order as a LinQuad
	// The opposite face should then have the remaining nodes in the same cyclic order
	// First find one face by computing max solid angle triplets
	std::vector<std::size_t> triplet0 =  Hexahedron::getNeighborNodes(0);
	assert(triplet0.size() == 3);
	std::vector<std::size_t> triplet1 =  Hexahedron::getNeighborNodes(triplet0[0]);
	std::vector<std::size_t> triplet2 =  Hexahedron::getNeighborNodes(triplet0[1]);
	assert(triplet1.size() == 3);
	assert(triplet2.size() == 3);
	// Now, find common node between triplet1 and triplet2
	std::size_t commonNodeID = triplet2.size();
	for(std::size_t k1 = 0; k1 < triplet1.size() && commonNodeID == triplet2.size(); ++k1)
	{
		for(std::size_t k2 = 0; k2 < triplet2.size() && commonNodeID == triplet2.size(); ++k2)
		{
			if((triplet1[k1] == triplet2[k2]) && (triplet1[k1] != 0))
				commonNodeID = triplet1[k1];
		}
	}
	std::vector<Point3D*> curPts(4);
	curPts[0] = nodeVec[0];
	curPts[1] = nodeVec[triplet0[0]];
	curPts[2] = nodeVec[commonNodeID];
	curPts[3] = nodeVec[triplet0[1]];
	// Get remaining 4 nodes, starting with remaining node from original triplet
	std::vector<Point3D*> curPts2(1, nodeVec[triplet0[2]]);
	curPts2.reserve(4);
	for(std::size_t k = 0; k < nodeVec.size(); ++k)
	{
		std::vector<Point3D*>::const_iterator fit1 = std::find(curPts.begin(), curPts.end(), nodeVec[k]);
		std::vector<Point3D*>::const_iterator fit2 = std::find(curPts2.begin(), curPts2.end(), nodeVec[k]);
		if(fit1 == curPts.end() && fit2 == curPts2.end())
			curPts2.push_back(nodeVec[k]);
	}
	assert(curPts2.size() == 4); // Fails if the 4 remaining points weren't found
	// Set ordering of second face, abusing LinearQuad to get correct ordering for 4 nodes
	LinearQuadrilateral<Point3D> tmpLinQuad(curPts2, Topologies::GenericMaterial());
	for(unsigned k = 0; k < curPts2.size(); ++k)
		curPts2[k] = tmpLinQuad.getNode(k);
	// Check that the 2 face normals point corrrectly
	Point3D v1 = *curPts[1] - *curPts[0], v2 = *curPts[2] - *curPts[0];
	Point3D cp1 = v1.crossProduct(v2);
	v1 = *curPts2[1] - *curPts2[0];
	v2 = *curPts2[2] - *curPts2[0];
	Point3D cp2 = v1.crossProduct(v2);
	if((cp1 * cp2) < 0) // Facing opposite directions
		std::reverse(curPts.begin() + 1, curPts.end());
	// Finally copy into one vector
	nodeVec.resize(curPts.size() + curPts2.size()); // Just in case
	for(std::size_t k = 0; k < curPts.size(); ++k)
		nodeVec[k] = curPts[k];
	for(std::size_t k = 0; k < curPts2.size(); ++k)
		nodeVec[k + curPts.size()] = curPts2[k];
	// Note that this does not check that the Jacobian is positive, that is done later!
}

Hexahedron::~Hexahedron()
{
}

Element<Point3D>* Hexahedron::getPatch(unsigned short int patchNum) const
{
	assert(patchNum < 6);
	return patchVec[patchNum];
}

Point3D* Hexahedron::getNode(unsigned short int nodeNum) const
{
	assert(nodeNum < 8);
	return nodeVec[nodeNum];
}

bool Hexahedron::isBoundaryCell() const
{
	Element<Point3D>* tempPatch;
	bool temp = false;
	for (unsigned int i = 0; i < 6; i++)
	{
		tempPatch = this->getPatch(i);
		temp |= tempPatch->isBoundaryPatch();
	}
	return temp;
}

std::unique_ptr<Element<Point3D>> Hexahedron::createPatch(unsigned short int patchNum) const
{
	assert(patchNum < 6);
	std::vector<Point3D*> facePtVec(4);
	switch(patchNum)
	{
	case 0:
		facePtVec[0] = nodeVec[0];
		facePtVec[1] = nodeVec[1];
		facePtVec[2] = nodeVec[2];
		facePtVec[3] = nodeVec[3];
		break;
	case 1:
		facePtVec[0] = nodeVec[0];
		facePtVec[1] = nodeVec[4];
		facePtVec[2] = nodeVec[5];
		facePtVec[3] = nodeVec[1];
		break;
	case 2:
		facePtVec[0] = nodeVec[0];
		facePtVec[1] = nodeVec[3];
		facePtVec[2] = nodeVec[7];
		facePtVec[3] = nodeVec[4];
		break;
	case 3:
		facePtVec[0] = nodeVec[6];
		facePtVec[1] = nodeVec[5];
		facePtVec[2] = nodeVec[4];
		facePtVec[3] = nodeVec[7];
		break;
	case 4:
		facePtVec[0] = nodeVec[6];
		facePtVec[1] = nodeVec[7];
		facePtVec[2] = nodeVec[3];
		facePtVec[3] = nodeVec[2];
		break;
	case 5:
		facePtVec[0] = nodeVec[6];
		facePtVec[1] = nodeVec[2];
		facePtVec[2] = nodeVec[1];
		facePtVec[3] = nodeVec[5];
		break;
	}
	return std::unique_ptr<Element<Point3D>>(new LinearQuadrilateral<Point3D>(facePtVec, itsMaterial));
}

void Hexahedron::addPatch(Element<Point3D>* thePatch)
{
	std::size_t freePatch = 0;
	while (patchVec[freePatch] != nullptr && freePatch < patchVec.size())
		++freePatch;
	assert(freePatch < patchVec.size());
	patchVec[freePatch] = thePatch;
}

void Hexahedron::switchXi2Xi3()
{
	std::reverse(nodeVec.begin() + 1, nodeVec.begin() + 4);
	std::reverse(nodeVec.begin() + 5, nodeVec.end());
}

void Hexahedron::setStandardPatchOrder()
{
	// Find patch that contains first 4 nodes
	std::vector<Point3D*> patchNodeVec(4);
	for(std::size_t k = 0; k < patchNodeVec.size(); ++k)
		patchNodeVec[k] = nodeVec[k];
	std::vector<Element<Point3D>*> newPatVec(6, nullptr);
	newPatVec[0] = findPatchByNodes(patchNodeVec);
	// Get patch that has last 4 nodes
	for(std::size_t k = 0; k < patchNodeVec.size(); ++k)
		patchNodeVec[k] = nodeVec[k + 4];
	newPatVec[1] = findPatchByNodes(patchNodeVec);
	// get 0, 1, 5, 4 patch (left side)
	patchNodeVec[0] = nodeVec[0];
	patchNodeVec[1] = nodeVec[1];
	patchNodeVec[2] = nodeVec[5];
	patchNodeVec[3] = nodeVec[4];
	newPatVec[2] = findPatchByNodes(patchNodeVec);
	// get 0, 3, 7, 4 patch (back side)
	patchNodeVec[0] = nodeVec[0];
	patchNodeVec[1] = nodeVec[3];
	patchNodeVec[2] = nodeVec[7];
	patchNodeVec[3] = nodeVec[4];
	newPatVec[3] = findPatchByNodes(patchNodeVec);
	// get 3, 7, 6, 2 patch (right side)
	patchNodeVec[0] = nodeVec[3];
	patchNodeVec[1] = nodeVec[7];
	patchNodeVec[2] = nodeVec[6];
	patchNodeVec[3] = nodeVec[2];
	newPatVec[4] = findPatchByNodes(patchNodeVec);
	// get 1, 2, 6, 5 patch (front side)
	patchNodeVec[0] = nodeVec[1];
	patchNodeVec[1] = nodeVec[2];
	patchNodeVec[2] = nodeVec[6];
	patchNodeVec[3] = nodeVec[5];
	newPatVec[5] = findPatchByNodes(patchNodeVec);
	if(std::find(newPatVec.begin(), newPatVec.end(), nullptr) == newPatVec.end())
		patchVec = newPatVec;
	else
	{
		std::cout << "Error in Hexahedron::setStandardPatchOrder, couldn't set order." << std::endl;
		std::cout << "Num unfound patches: " << std::count(newPatVec.begin(), newPatVec.end(), nullptr) << std::endl;
		std::cout << "Num uncreated patches: " << std::count(patchNodeVec.begin(), patchNodeVec.end(), nullptr) << std::endl;
		std::cout << "pts = [";
		for(unsigned k = 0; k < 8; ++k)
			std::cout << nodeVec[k] << ": " << *nodeVec[k] << std::endl;;
		std::cout << "];" << std::endl;
		std::cout << "Patches: " << std::endl;
		for(unsigned k = 0; k < 6; ++k)
		{
			std::cout << k << ": pts: " << std::endl;
			for(unsigned k2 = 0; k2 < patchVec[k]->getNumNodes(); ++k2)
				std::cout << patchVec[k]->getNode(k2) << ": " << *patchVec[k]->getNode(k2) << std::endl;
		}
		std::vector<Point3D*> tstPts(4);
		for(unsigned k = 0; k < 4; ++k)
			tstPts[k] = nodeVec[k];
		std::cout << "TESTING LINEARQUAD:::::" << std::endl;
		LinearQuadrilateral<Point3D> tmpLinQuad(tstPts, Topologies::GenericMaterial(), true);
		std::cout << "volume: " << volumeIntegral() << std::endl;
	}
}

Element<Point3D>* Hexahedron::findPatchByNodes(const std::vector<Point3D*>& patchNodeVec) const
{
	for(std::size_t k = 0; k < patchVec.size(); ++k)
	{
		std::vector<Point3D*> curPatchNodeVec = getPatchNodeVec(k);
		if(std::is_permutation(patchNodeVec.begin(), patchNodeVec.end(), curPatchNodeVec.begin()))
			return patchVec[k];
	}
	return nullptr;
}

std::vector<Point3D*> Hexahedron::getPatchNodeVec(std::size_t kp) const
{
	std::vector<Point3D*> patchNodeVec(patchVec[kp]->getNumNodes());
	for(std::size_t kn = 0; kn < patchNodeVec.size(); ++kn)
		patchNodeVec[kn] = patchVec[kp]->getNode(kn);
	return patchNodeVec;
}

void Hexahedron::addPatchBases(Element<Point3D>* theCPatch, Point3D interpPos[], std::size_t startPoint, int signum)
{
	// Not implemented, needs to be for higher order elements
}

bool Hexahedron::containsNode(Point3D* tstNode) const
{
	std::vector<Point3D*>::const_iterator fit = std::find(nodeVec.begin(), nodeVec.end(), tstNode);
	return fit != nodeVec.end();
}

void Hexahedron::addNodeBases(const Point3D* theNode, std::size_t bfGID)
{
	bool found = false;
	for(std::size_t k = 0; k < nodeVec.size() && !found; ++k)
	{
		if(theNode == nodeVec[k])
		{
			bfVec[k] = bfGID;
			found = true;
		}
	}
	assert(found); // Fails if theNode was not found
	for(unsigned k = 0; k < patchVec.size(); k++)
	{
		if(patchVec[k]->containsNode(theNode))
			patchVec[k]->addNodeBases(theNode, bfGID);
	}
}

void Hexahedron::addCellBases(std::size_t& lastAdded)
{
	// Not implemented, should be for higher order elements
}

