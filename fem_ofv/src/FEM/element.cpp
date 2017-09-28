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

#include "element.h"
#include "elemedge.h"
#include <algorithm>
#include "point3d.h"

template <typename T>
Element<T>::Element(const std::vector<T*>& ptVec, const Topologies::GenericMaterial& inMat):
	nodeVec(ptVec),
	material(inMat),
	reordered(false),
	pCell0(nullptr),
	pCell1(nullptr)
{
	edgeVec.reserve(nodeVec.size());
	for(unsigned k = 0; k < nodeVec.size(); ++k)
		edgeVec.push_back(nullptr);
	bfVec.resize(nodeVec.size());
	for(unsigned short int i = 0; i < bfVec.size(); i++)
		bfVec[i] = 0;
}

template <typename T>
Element<T>::~Element()
{
}

template <typename T>
bool Element<T>::operator==(const Element<T>& arg2) const
{
	if(myPT != arg2.myPT)
		return false;
	if(nodeVec.size() != arg2.nodeVec.size())
		return false;
	return std::is_permutation(nodeVec.begin(), nodeVec.end(), arg2.nodeVec.begin());
}

template <typename T>
bool Element<T>::addCell(Cell* inCell) 
{
	if (pCell0 == nullptr)
		pCell0 = inCell;
	else
	{
		if(pCell1 != nullptr) // Check if element already has 2 associated Cells
			return true;
		pCell1 = inCell;
	}
	return false;
}

template <typename T>
std::unique_ptr<ElemEdge<T>> Element<T>::createEdge(unsigned short int edgeNum) const
{
	assert(edgeNum < nodeVec.size());
	return std::unique_ptr<ElemEdge<T>>(new ElemEdge<T>(nodeVec[edgeNum], nodeVec[(edgeNum + 1) % nodeVec.size()]));
}

template <typename T>
void Element<T>::addEdge(ElemEdge<T>* theEdge)
{
	unsigned short int freeEdge = 0;
	while (edgeVec[freeEdge] != nullptr && freeEdge < edgeVec.size())
		++freeEdge;
	assert(freeEdge < edgeVec.size());
	edgeVec[freeEdge] = theEdge;
}

template <typename T>
void Element<T>::switchXiEta()
{
	std::reverse(std::next(nodeVec.begin()),nodeVec.end());
}

template <typename T>
void Element<T>::setStandardEdgeOrder()
{
	if(edgeVec.size() < 3)
		return;
	// Note that this assumes that the nodes are in cyclic order
	std::vector<ElemEdge<T>*> tempEdge(edgeVec.size(), nullptr);
	for(std::size_t iNode = 0; iNode < edgeVec.size(); ++iNode)
	{
		T* pNode0 = nodeVec[iNode];
		T* pNode1 = nodeVec[(iNode + 1)%edgeVec.size()];
		for(std::size_t jEdge = 0; jEdge < edgeVec.size(); ++jEdge)
		{
			if(!(edgeVec[jEdge]->containsNode(pNode0)) && edgeVec[jEdge]->containsNode(pNode1))
				tempEdge[iNode] = edgeVec[jEdge];
		}
	}
	// Check for nullptr and copy
	if(std::find(tempEdge.begin(), tempEdge.end(), nullptr) == tempEdge.end())
		edgeVec = tempEdge;
	else
		std::cout << "Warning, setStandardEdgeOrder didn't work!  Nodes may be unordered" << std::endl;
}

template <typename T>
bool Element<T>::containsNode(const T* tstNode) const
{
	typename std::vector<T*>::const_iterator fit = std::find(nodeVec.begin(), nodeVec.end(), tstNode);
	return fit != nodeVec.end();
}

template <typename T>
void Element<T>::addNodeBases(const T* theNode, std::size_t bfGID)
{
	bool found = false;
	for(unsigned k = 0; k < nodeVec.size() && !found; ++k)
	{
		if(theNode == nodeVec[k])
		{
			bfVec[k] = bfGID;
			found = true;
		}
	}
	assert(found);
	for(unsigned k = 0; k < edgeVec.size(); k++)
	{
		if(edgeVec[k])
		{
			if(edgeVec[k]->containsNode(theNode))
				edgeVec[k]->addNodeBases(theNode, bfGID);
		}
	}
}

template <typename T>
std::size_t Element<T>::getGlobalBFonBoundary(unsigned locBF) const
{
	ElemEdge<T>* bdEdge = findBoundaryEdge();
	return bdEdge->getGlobalBF(locBF);
}

template <typename T>
ElemEdge<T>* Element<T>::findBoundaryEdge() const
{
	ElemEdge<T>* bdryEdge = nullptr;
	for(unsigned k = 0; k < edgeVec.size() && bdryEdge == nullptr; ++k)
		if(edgeVec[k]->isBoundaryEdge())
			bdryEdge = edgeVec[k];
	assert(bdryEdge != nullptr);
	return bdryEdge;
}

template <typename T>
T Element<T>::centroid() const
{
	T tmp; // Point2D and Point3D initialize all components to 0.
	for(unsigned k = 0; k < nodeVec.size(); ++k)
		tmp += *(nodeVec[k]);
	tmp /= (double)nodeVec.size();
	return tmp;
}

template <typename T>
unsigned Element<T>::getLocEdgeID(ElemEdge<T>* findEdge) const
{
	unsigned id = 0;
	for(unsigned k = 0; k < edgeVec.size(); k++)
	{
		if(findEdge == edgeVec[k])
			id = k;
	}
	return id;
}

template <typename T>
std::pair<double,double> Element<T>::getCijs(const ElasticProblemType inEPT) const
{
	double rho = this->material.getParameter(cRho), lambda = this->material.getParameter(cLambda),
			 mu = this->material.getParameter(cMu);
	double cii, cij;
	if(inEPT == eptPlaneStrain)
	{
		cii = lambda + 2.*mu;
		cij = lambda;
	}
	else
	{
		cii = 4.*mu*(lambda + mu)/(lambda + 2.*mu);
		cij = 2*lambda*mu/(lambda + 2.*mu);
	}
	return std::make_pair(cii, cij);
}

template class Element<Point2D>;
template class Element<Point3D>;
