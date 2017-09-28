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

#ifndef CSGNODE_H
#define CSGNODE_H

#include <vector>
#include <iostream>
#include "cgal_types.h"
#include "genericmaterial.h"
#include "inputloader.h"

namespace Topologies{
//! An enumeration to define the type of material to use
enum DecodeType { dtBoolean, dtRandMat, dtDatabaseMat };
//! An enumeration to define the type of geometric function, union or subtraction
enum FunctionType {ftUnion = 0, ftSubtraction};
//! An enuperation to define which coordinates of a point can be modified
enum FixPointCoordType { fpctXYZ, fpctXY, fpctX, fpctY, fpctXZ, fpctYZ, fpctZ};
//! An enumeration to define different types of mirroring
enum MirrorType{mtNone, mtX, mtY, mtZ};

//! An abstract class that defines constructive solid geometry (CSG) nodes
class CSGNode
{
public:
	//! Constructor that takes the current depth in the tree
	CSGNode(int curTreeDepth);
	virtual ~CSGNode();

	//! Generates a 2d topology by evaluating the subtree rooted at this CSGNode
	/*! This function will be called recursively to generate the entire topology */
	virtual void decode(Nef_polyhedron_2& inNP) const = 0;
	//! Generates a 2d topology mask by evaluating the subtree rooted at this CSGNode
	/*! The mask gives the topology defined by computing the union of all separate regions,
	 *   i.e. those with different materials.
	 */
	virtual void decodeMask(Nef_polyhedron_2& inNP) const = 0;
	//! Determines the material properties at the point `testPt`
	virtual bool decodeMaterial(Point_2& testPt, GenericMaterial& theGM) const = 0;
	//! Generates a 3d topology by evaluating the subtree rooted at this CSGNode
	/*! This function will be called recursively to generate the entire topology */
	virtual void decode(Nef_polyhedron_3& inNP) const = 0;
	//! Generates a 3d topology mask by evaluating the subtree rooted at this CSGNode
	/*! The mask gives the topology defined by computing the union of all separate regions,
	 *   i.e. those with different materials.
	 */
	virtual void decodeMask(Nef_polyhedron_3& inNP) const = 0;
	//! Determines the material properties at the point `testPt`
	virtual bool decodeMaterial(Point_3& testPt, GenericMaterial& theGM) const = 0;
	//! Determines the ID of the material at the point `testPt`
	virtual bool decodeMaterial(Point_3& testPt, int& theGM) const = 0;
	//! Generates a vector of GenericMaterial objects that contains all materials in the subtree rooted here
	virtual void setMaterialList(std::vector<GenericMaterial>& theGMVec) const = 0;
	//! Returns the number of points contained in the subtree rooted at this node.  
	virtual unsigned countPoints() const = 0;
	//! Returns the number of tree nodes contained in the subtree rooted at this node.  
	virtual unsigned countNodes() const = 0;
	//! Writes the data of the subtree rooted here to the argument ofstream
	virtual void printCSGNodeToFile(std::ofstream& popFile) const = 0;
	//! Generates a set of vectors containing information about the subtree needed to reconstruct this subtree in 2d
	virtual void getMPIDataFormat(std::vector<int>& treeVec, std::vector<W_Point_2>& ptVec, 
				 std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const = 0;
	//! Generates a set of vectors containing information about the subtree needed to reconstruct this subtree in 3d
	virtual void getMPIDataFormat(std::vector<int>& treeVec, std::vector<Point_3>& ptVec,
                 std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const = 0;
	//! Determines whether or not a subtree is safe to mesh by checking some geometric features
	/*! Testing of the CGAL mesher has shown it to crash in some circumstances, such as when points are too
	 *  close together.  In addition, points and segments that are too close generate large numbers of elements.
	 */
	virtual bool checkMeshability() = 0;
	//! Returns the maximum radius of the points in the subtree rooted here
	virtual Real getMaxPointRadius() const = 0;
	//! Translates all points contained in the subtree rooted here outside of the bounds set in the input file
	virtual void movePointsOB() = 0;
	//! Sets the coordinates of all points contained in this subtree to the values in the iterator arguments
	virtual void setVector(std::vector<double>::const_iterator& curPos, const std::vector<double>::const_iterator& endPos) = 0;
	//! Adds all points in the subtree rooted here to the list `ptList`
	virtual void getLocalOptFormat(std::list<W_Point_2>& ptList) const = 0;
	//! Returns the number of children in the subtree rooted here
	virtual unsigned getNumChildren() const = 0;
	//! Genertes a vector of the CSGTerminalNode (or CSGTerminalNode3D) objects in this subtree
	virtual void getTerminalNodeList(std::vector<CSGNode*>& nodeVec) = 0;

	//! Sets the tree depth (number of parents) of this node to `inDepth`
	void setDepthInTree(int inDepth);
	//! Returns the parent of this node
	CSGNode* getParent();
	//! Returns the tree depth (number of parents) of this node
	int getDepthInTree();

	//! Set the DecodeType, or how materials are handled for all CSGTree
	static void setDecodeType(DecodeType inDT);
	//! Sets the DimensionType (2 or 3)
	static void setDimensionType(DimensionType inDT);
	//! Sets a value within which points are snapped directly to the boundary
	static void setBoundarySnap(Real inBoundSnap);
	//! Returns the current DecodeType
	static DecodeType getDecodeType();
	//! Returns the current DimensionType
	static DimensionType getDimensionType();
	//! Returns the list of GenericMaterial objects
	static std::vector<GenericMaterial> getMaterialList();
	//! Returns the current boundary snap value
	static Real getBoundarySnap();
protected:
	CSGNode* parent;
	int depthInTree;
	static DecodeType theDT;
	static DimensionType theDim;
	static std::vector<GenericMaterial> materialList;
	static int numFuncTypes;
	static Real boundMag, boundSnap, prioritymax;
	static bool useBoundSnap;
protected:
	CSGNode(const CSGNode& copy);
	CSGNode(CSGNode && copy);
	CSGNode& operator=(const CSGNode& copy);
	void swap(CSGNode& arg2);
};

inline
CSGNode* CSGNode::getParent()
{
	return parent;
}

inline
int CSGNode::getDepthInTree()
{
	return depthInTree;
}

inline
void CSGNode::setDepthInTree(int inDepth)
{
	depthInTree = inDepth;
}

inline
void CSGNode::setDecodeType(DecodeType inDT)
{
	theDT = inDT;
	numFuncTypes = 2;
}

inline
void CSGNode::setDimensionType(DimensionType inDT)
{
    theDim = inDT;
}

inline
DecodeType CSGNode::getDecodeType()
{
	return theDT;
}

inline
DimensionType CSGNode::getDimensionType()
{
    return theDim;
}

inline
std::vector<GenericMaterial> CSGNode::getMaterialList()
{
	return materialList;
}

inline
void CSGNode::setBoundarySnap(Real inBoundSnap)
{
	boundSnap = inBoundSnap;
	useBoundSnap = true;
}

inline 
Real CSGNode::getBoundarySnap()
{
	return boundSnap;
}
} //namespace

#endif

