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

#ifndef CSGFUNCTIONNODE_H
#define CSGFUNCTIONNODE_H

#include <vector>
#include <iostream>
#include "csgnode.h"

namespace Topologies{
//! A class that defines function nodes (non-terminal nodes) for a constructive solid geometry tree
/*! This implementation of CSGNode defines function nodes, those that perform geometric Boolean
 *  operations on CSGTerminalNode objects.  Two functions are defined, union and subtraction
 */
class CSGFunctionNode : public CSGNode
{
public:
	//! Constructor that takes the number of CSGNode objects to generate along 2 axes, along with the number of points and the regions height and width
	CSGFunctionNode(int curTreeDepth, unsigned numNodesX, unsigned numNodesY, unsigned numPts, double width, double height);
	CSGFunctionNode(int curTreeDepth, const std::vector<std::list<W_Point_2>>& ptListVec);
	CSGFunctionNode(const CSGFunctionNode& copy);
	CSGFunctionNode(CSGFunctionNode && copy);
	CSGFunctionNode& operator=(CSGFunctionNode copy);
	void swap(CSGFunctionNode& arg2);
	//! A specialized copy constructor that uses a MirrorType to reflect the subtree's nodes in a specified way
	CSGFunctionNode(const CSGFunctionNode& copy, MirrorType inMT);
	//! Constructor that loads a CSGNode from a file
	CSGFunctionNode(int curTreeDepth, int functionType, std::ifstream& saveFile);
	//! Constructor that sets a subtree up with the specified values in C-style arrays
	/*! 
	 * @param curTreeDepth Current depth in tree
	 * @param numMatParams Number of material properties 
	 * @param functionType The FunctionType of the node, cast as an int
	 * @param ptVecPos Current position in the `ptArray` array
	 * @param matPos Current position in the `matArray` array
	 * @param prPos Current position in the `prArray`
	 * @param treeVecPos Current position in the `treeVecPos`
	 * @param ptArray Array of coordinates and weighted alpha shape values
	 * @param matArray Array of material properties
	 * @param prArray Array of alpha values for computing alpha shapes of terminal nodes
	 * @param treeArray Array containing information about tree such as function node types and numbers of points in terminal nodes
	 */
	CSGFunctionNode(int curTreeDepth, int numMatParams, int functionType, int& ptVecPos, int& matPos, int& prPos, 
				int& treeVecPos, Real* ptArray, Real* matArray, Real* prArray, int* treeArray);
	//! Constructor that sets a subtree up with the specified values, set in vector iterators
	/*! An implementation of the C-style constructor above, uses iterators */
	CSGFunctionNode(int curTreeDepth, int functionType, std::vector<int>::const_iterator& nmIt, 
									std::vector<double>::const_iterator& ptIt, std::vector<double>::const_iterator& matIt, 
									std::vector<double>::const_iterator& prIt, std::vector<int>::const_iterator& treeIt);
	//! Constructor that takes ownership the children nodes in `inKids`
	CSGFunctionNode(int curTreeDepth, FunctionType functionType, CSGNode* parent, std::vector<std::unique_ptr<CSGNode> >& inKids);
	virtual ~CSGFunctionNode();

	virtual void decode(Nef_polyhedron_2& inNP) const;
	virtual void decodeMask(Nef_polyhedron_2& inNP) const;
	virtual bool decodeMaterial(Point_2& testPt, GenericMaterial& theGM) const;
	virtual void decode(Nef_polyhedron_3& inNP) const;
	virtual void decodeMask(Nef_polyhedron_3& inNP) const;
	virtual bool decodeMaterial(Point_3& testPt, GenericMaterial& theGM) const;
	virtual bool decodeMaterial(Point_3& testPt, int& theGM) const;
	virtual void setMaterialList(std::vector<GenericMaterial>& theGMVec) const;
	virtual unsigned countPoints() const;
	virtual unsigned countNodes() const;
	virtual void printCSGNodeToFile(std::ofstream& popFile) const;
	virtual void getMPIDataFormat(std::vector<int>& treeVec, std::vector<W_Point_2>& ptVec, 
				  std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const;
	virtual void getMPIDataFormat(std::vector<int>& treeVec, std::vector<Point_3>& ptVec,
                  std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const;
	virtual Real getMaxPointRadius() const;
	virtual bool checkMeshability();
	virtual void movePointsOB();
	virtual void getLocalOptFormat(std::list<W_Point_2>& ptList) const;
	virtual void setVector(std::vector<double>::const_iterator& curPos, const std::vector<double>::const_iterator& endPos);
	virtual unsigned getNumChildren() const;
	virtual void getTerminalNodeList(std::vector<CSGNode*>& nodeVec);

	//! Returns the FunctionType of this CSGFunctionNode
	FunctionType getFunctionType();
	//! Sets the FunctionType of this CSGFunctionNode
	void setFunctionType(FunctionType inFT);

private:
	void createChildren(unsigned numNodes);
	void createChildren(unsigned numNodesX, unsigned numNodesY, unsigned numPts, double width, double height);
	void unionFunction(std::vector<Nef_polyhedron_2>& sdvec, Nef_polyhedron_2& inNP) const;
	void subtractionFunction(std::vector<Nef_polyhedron_2>& sdvec, Nef_polyhedron_2& inNP) const;
	void subtractionFunction(std::vector<Nef_polyhedron_2>& sdvec, Nef_polyhedron_2& rop, Nef_polyhedron_2& inNP) const;
	void overlapFunction(std::vector<Nef_polyhedron_2>& sdvec_exc, Nef_polyhedron_2& rop, Nef_polyhedron_2& inNP) const;
	void unionFunction(std::vector<Nef_polyhedron_3>& sdvec, Nef_polyhedron_3& inNP) const;
	void subtractionFunction(std::vector<Nef_polyhedron_3>& sdvec, Nef_polyhedron_3& inNP) const;
	void subtractionFunction(std::vector<Nef_polyhedron_3>& sdvec, Nef_polyhedron_3& rop, Nef_polyhedron_3& inNP) const;
	void overlapFunction(std::vector<Nef_polyhedron_3>& sdvec_exc, Nef_polyhedron_3& rop, Nef_polyhedron_3& inNP) const;

private:
	std::vector<std::unique_ptr<CSGNode> > children;
	FunctionType theFT;
};

inline 
void CSGFunctionNode::getTerminalNodeList(std::vector<CSGNode*>& nodeVec)
{
	for(unsigned k = 0; k < children.size(); ++k)
		children[k]->getTerminalNodeList(nodeVec);
}

inline
unsigned CSGFunctionNode::getNumChildren() const
{
  return children.size();
}

inline
void CSGFunctionNode::setVector(std::vector<double>::const_iterator& curPos, const std::vector<double>::const_iterator& endPos)
{
	for(Uint k = 0; k < children.size(); ++k)
		children[k]->setVector(curPos, endPos);
}

inline
void CSGFunctionNode::getLocalOptFormat(std::list<W_Point_2>& ptList) const
{
	for(Uint k = 0; k < children.size(); ++k)
		children[k]->getLocalOptFormat(ptList);
}

inline
FunctionType CSGFunctionNode::getFunctionType()
{
	return theFT;
}

inline
void CSGFunctionNode::setFunctionType(FunctionType inFT)
{
	theFT = inFT;
}
} // namespace
#endif

