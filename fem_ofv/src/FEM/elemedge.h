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

#ifndef TRIEDGE_H
#define TRIEDGE_H

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <functional>
#include "point2d.h"

template <typename T> class Element;

//! Class implements an edge adjoining two Element objects and containts two defining nodes
/*! This class is used for connectivity (topological) information about a mesh.  It is used
*  to define a unique edge between two Element objects.  The template parameter gives the
*  type of point to use for the nodes, Point2D or Point3D.
*/
template <class T>
class ElemEdge 
{
public:
	//! Default constructor initializes all pointers to nullptr
	ElemEdge();
	//! Constructor sets the two defining nodes of the ElemEdge
	ElemEdge(T* pInNode1, T* pInNode2);
	ElemEdge<T>& operator=(ElemEdge<T> rhs);
	void swap(ElemEdge<T>& edge2);
	~ElemEdge();

	//! Equality Test, edges with any permutation of nodes are considered equal
	bool operator==(const ElemEdge<T>& rhs) const;
	//! Adds a new element to the edge.
	/*! The first time this function is called, the argument is stored in pElem0 and the second time it is stored in pElem1.
	 *  Two elements are sufficient for 2d.  In 3d, subsequent elements are stored in elemVec.
	 */
	bool addElement(Element<T>* theElem);
	//! Returns the first element associated with this ElemEdge
	Element<T>* elem0() const;
	//! Returns the second element associated with this ElemEdge
	Element<T>* elem1() const;
	//! Returns the Nth element associated with this ElemEdge.
	/*! The indexing is such that N=0 gives elem0, N=1 gives elem1 and N>1 gives elemVec[N-2]*/
	Element<T>* elemN(unsigned N) const;
	//! Returns the number of Elements this ElemEdge is associated with
	/*! Note that this is only valid (and used) in 3d. */
	unsigned getNumAttachedElems() const {return 2 + elemVec.size();}
	//! Returns the first node
	T* node0() const;
	//! Returns the second node
	T* node1() const;
	//! Reterns whether or not this ElemEdge is a boundary edge.  Only valid in 2d.
	/*! The test is whether or not there is only 1 element associated with this edge. */
	bool isBoundaryEdge() const;
	//! Returns whether or not this ElemEdge is equal to theEdge, but has its nodes in opposite order.
	bool isAlignedWith(const ElemEdge<T>& theEdge) const;
	//! Returns whether or not theNode is contained in this ElemEdge
	bool containsNode(const T* theNode) const;
	//! Returns the physical length of the edge
	double length() const;
	//! Returns a vector from one point to the other
	T edgeVec() const;
	//! Adds nodal basis function id to this ElemEdge
	void addNodeBases(const T* theNode, std::size_t bfGID);
	//! Returns the global basis function id associated with local id locBF
	std::size_t getGlobalBF(unsigned locBF) const;
private:
	Element<T>* pElem0, *pElem1;
	std::vector<Element<T>*> elemVec;
	T* pNode0, *pNode1;
	std::vector<std::size_t> gbfs;
};

//! Hash function for testing ElemEdge equality, used in std::unordered_set
/*! Adapted from http://stackoverflow.com/questions/1536393/good-hash-function-for-permutations 
*  For this class, we want to ensure that edges that contain permutations of nodes test
*  as equal, so a hash that gives the same value for any permutation of node pointers is needed.
*/
template<typename T>
struct Elemedge_hash
{
	typedef const ElemEdge<T>* const argument_type;
	typedef std::size_t result_type;
	result_type operator()(argument_type const& s) const
	{
		const std::size_t R = 1779033703;
		std::size_t h1 = reinterpret_cast<size_t>(s->node0());
		std::size_t h2 = reinterpret_cast<size_t>(s->node1());
		return (R + 2*h1)*(R + 2*h2)/2;
	}
};

//! ElemEdge equality testing struct, used in std::unordered_set
template<typename T>
struct Elemedge_equal
{
	typedef const ElemEdge<T>* const argument_type;
	bool operator()(argument_type const& x, argument_type const& y) const
	{
		return *x == *y;
	}
};

template <typename T> 
Element<T>* ElemEdge<T>::elem0() const
{
	return pElem0;
}

template <typename T> 
Element<T>* ElemEdge<T>::elem1() const 
{
	return pElem1;
}

template <typename T> 
Element<T>* ElemEdge<T>::elemN(unsigned N) const
{
	if(N == 0)
		return pElem0;
	else if(N == 1)
		return pElem1;
	N -= 2;
	assert(N < elemVec.size());
	return elemVec[N];
}

template <typename T>
T* ElemEdge<T>::node0() const
{
	return pNode0;
}

template <typename T>
T* ElemEdge<T>::node1() const 
{
	return pNode1;
}

template <typename T> 
bool ElemEdge<T>::isBoundaryEdge() const
{
	return pElem1 == nullptr;
}

template <typename T>
bool ElemEdge<T>::isAlignedWith(const ElemEdge<T>& theEdge) const
{
	return (pNode1 == theEdge.pNode0) && (pNode0 == theEdge.pNode1);
}

template <typename T>
bool ElemEdge<T>::containsNode(const T* theNode) const
{
	return (pNode0 == theNode) || (pNode1 == theNode);
}

template <typename T>
double ElemEdge<T>::length() const
{
	return abs(*pNode0 - *pNode1);
}

template <typename T>
T ElemEdge<T>::edgeVec() const
{
	T tmpvec(*pNode0 - *pNode1);
	return tmpvec;
}

template <typename T>
std::size_t ElemEdge<T>::getGlobalBF(unsigned locBF) const
{
	assert(locBF < 2);
	return gbfs[locBF];
}

#endif

