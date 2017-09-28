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

#include "tetrahedron.h"

Tetrahedron::Tetrahedron(const std::vector<Point3D*>& inNodeVec, const Topologies::GenericMaterial& inMat):
	Cell(inMat)
{
	assert(inNodeVec.size() == 4);
	nodeVec = inNodeVec;
	patchVec = std::vector<Element<Point3D>*>(4, nullptr);
	bfVec.resize(nodeVec.size());
	for(std::size_t i = 0; i < bfVec.size(); i++)
		bfVec[i] = 0;
}

Tetrahedron::~Tetrahedron()
{
}

Element<Point3D>* Tetrahedron::getPatch(unsigned short int patchNum) const
{
	assert(patchNum < 4);
	return patchVec[patchNum];
}

Point3D* Tetrahedron::getNode(unsigned short int nodeNum) const
{
	assert(nodeNum < 4);
	return nodeVec[nodeNum];
}

bool Tetrahedron::isBoundaryCell() const
{
	Element<Point3D>* tempPatch;
	bool temp = false;
	for (unsigned int i = 0; i < 4; i++)
	{
		tempPatch = this->getPatch(i);
		temp |= tempPatch->isBoundaryPatch();
	}
	return temp;
}

std::unique_ptr<Element<Point3D>> Tetrahedron::createPatch(unsigned short int patchNum) const
{
	assert(patchNum < 4);
	std::vector<Point3D*> facePtVec(3);
	facePtVec[0] = nodeVec[patchNum];
	facePtVec[1] = nodeVec[(patchNum + 1) % 4];
	facePtVec[2] = nodeVec[(patchNum + 2) % 4];
	return std::unique_ptr<Element<Point3D>>(new LinearTriangle<Point3D>(facePtVec, itsMaterial));
}

void Tetrahedron::addPatch(Element<Point3D>* thePatch)
{
	std::size_t freePatch = 0;
	while (patchVec[freePatch] != nullptr && freePatch < patchVec.size())
		++freePatch;
	assert(freePatch < patchVec.size());
	patchVec[freePatch] = thePatch;
}

void Tetrahedron::switchXi2Xi3()
{
	std::swap(nodeVec[2], nodeVec[3]);
}

void Tetrahedron::setStandardPatchOrder()
{
	Element<Point3D>* tempPatch[4];
	Point3D* theNode;
	for (unsigned iNode = 0; iNode < 4; iNode++)
	{
		theNode = nodeVec[iNode];
		for (unsigned jPatch = 0; jPatch < 4; jPatch++)
		{
			if (!(patchVec[jPatch]->containsNode(theNode)))
				tempPatch[iNode] = patchVec[jPatch];
		}
	}
	for (unsigned jPatch = 0; jPatch < 4; jPatch++)
		patchVec[jPatch] = tempPatch[jPatch];
}

void Tetrahedron::addPatchBases(Element<Point3D>* theCPatch, Point3D interpPos[], std::size_t startPoint, int signum)
{
	// Not implemented, needs to be for higher order elements
}

bool Tetrahedron::containsNode(Point3D* tstNode) const
{
	std::vector<Point3D*>::const_iterator fit = std::find(nodeVec.begin(), nodeVec.end(), tstNode);
	return fit != nodeVec.end();
}

void Tetrahedron::addNodeBases(const Point3D* theNode, std::size_t bfGID)
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
	assert(found); // Fails if theNode doesn't belong to this element
	for(std::size_t k = 0; k < patchVec.size(); k++)
	{
		if(patchVec[k]->containsNode(theNode))
			patchVec[k]->addNodeBases(theNode, bfGID);
	}
}

void Tetrahedron::addCellBases(std::size_t& lastAdded)
{
	// Not implemented, should be for higher order elements
}

