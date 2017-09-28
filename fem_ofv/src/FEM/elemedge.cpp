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

#include "elemedge.h"
#include "point3d.h"
#include "point2d.h"

template <typename T>
ElemEdge<T>::ElemEdge():
	pElem0(nullptr),
	pElem1(nullptr),
	pNode0(nullptr),
	pNode1(nullptr)
{
	gbfs.resize(2);
}

template <typename T>
ElemEdge<T>::ElemEdge(T* pInNode0, T* pInNode1):
	pElem0(nullptr),
	pElem1(nullptr),
	pNode0(pInNode0),
	pNode1(pInNode1)
{
	gbfs.resize(2);
}

template <typename T>
ElemEdge<T>::~ElemEdge()
{
}

template <typename T>
ElemEdge<T>& ElemEdge<T>::operator=(ElemEdge<T> rhs)
{
	swap(rhs);
	return *this;
}

template <typename T>
bool ElemEdge<T>::operator==(const ElemEdge<T>& rhs) const
{
	bool temp1 = (pNode0 == rhs.pNode0) && (pNode1 == rhs.pNode1);
	bool temp2 = (pNode1 == rhs.pNode0) && (pNode0 == rhs.pNode1);
	return temp1 || temp2;
}

template <typename T>
bool ElemEdge<T>::addElement(Element<T>* theTri)
{
	if (pElem0 == 0)
		pElem0 = theTri;
	else
	{
		if(pElem1 != 0)
		{
			elemVec.push_back(theTri);
			return true;
		}
		pElem1 = theTri;
	}
	return false;
}

template <typename T>
void ElemEdge<T>::swap(ElemEdge<T>& edge2)
{
	std::swap(pElem0, edge2.pElem0);
	std::swap(pElem1, edge2.pElem1);
	std::swap(pNode0, edge2.pNode0);
	std::swap(pNode1, edge2.pNode1);
	elemVec.swap(edge2.elemVec);
	gbfs.swap(edge2.gbfs);
}

template <typename T>
void ElemEdge<T>::addNodeBases(const T* theNode, std::size_t bfGID)
{
	if(theNode == pNode0)
		gbfs[0] = bfGID;
	else if(theNode == pNode1)
		gbfs[1] = bfGID;
	else
	{
		std::cerr << "Node not found in edge for basis ID allocation." << std::endl;
		assert(false);
	}
}

template class ElemEdge<Point2D>;
template class ElemEdge<Point3D>;

