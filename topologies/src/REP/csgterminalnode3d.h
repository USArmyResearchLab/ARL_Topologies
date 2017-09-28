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

#ifndef CSGTERMINALNODE3D_H
#define CSGTERMINALNODE3D_H

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include "csgnode.h"
#include "genericmaterial.h"
#include "helper.h"

namespace Topologies{
class CSGTerminalNode3D : public CSGNode
{
public:
	CSGTerminalNode3D(int curTreeDepth);
	CSGTerminalNode3D(const CSGTerminalNode3D& copy);
	//! A specialized copy constructor that uses a MirrorType to reflect the subtree's nodes in a specified way
	CSGTerminalNode3D(const CSGTerminalNode3D& copy, MirrorType inMT);
	//! Constructor that loads a CSGNode from a file
	CSGTerminalNode3D(int curTreeDepth, std::ifstream& saveFile);
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
	CSGTerminalNode3D(int curTreeDepth, int numMatParams, int& ptVecPos, int& matPos, int& prPos, int& treeVecPos, 
			  Real* ptArray, Real* matArray, Real* prArray, int* treeArray);
	//! Constructor that sets a subtree up with the specified values, set in vector iterators
	/*! An implementation of the C-style constructor above, uses iterators */
	CSGTerminalNode3D(int curTreeDepth, int numMatParams, std::vector<double>::const_iterator& ptIt,
                    std::vector<double>::const_iterator& matIt, std::vector<double>::const_iterator& prIt,
                    std::vector<int>::const_iterator& treeIt);
	//! Constructor that copies the points in `ptList`
	CSGTerminalNode3D(std::list<Point_3>& ptList);
	CSGTerminalNode3D(CSGTerminalNode3D&& copy);
	CSGTerminalNode3D& operator=(CSGTerminalNode3D rhs);
	void swap(CSGTerminalNode3D& arg2); 
	virtual ~CSGTerminalNode3D();

	virtual void decode(Nef_polyhedron_3& inNP) const;
	virtual void decodeMask(Nef_polyhedron_3& inNP) const;
	virtual bool decodeMaterial(Point_3& testPt, GenericMaterial& theGM) const;
	virtual bool decodeMaterial(Point_3& testPt, int& theGM) const;
	virtual void decode(Nef_polyhedron_2& inNP) const;
	virtual void decodeMask(Nef_polyhedron_2& inNP) const;
	virtual bool decodeMaterial(Point_2& testPt, GenericMaterial& theGM) const;
	virtual void setMaterialList(std::vector<GenericMaterial>& theGMVec) const;
	virtual unsigned countPoints() const;
	virtual unsigned countNodes() const;
	virtual void printCSGNodeToFile(std::ofstream& popFile) const;
	virtual void getMPIDataFormat(std::vector<int>& treeVec, std::vector<Point_3>& ptVec, 
					std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const;
	virtual void getMPIDataFormat(std::vector<int>& treeVec, std::vector<W_Point_2>& ptVec,
                  std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const;
	virtual Real getMaxPointRadius() const;
	virtual bool checkMeshability();
	virtual void movePointsOB();
	virtual void setVector(std::vector<double>::const_iterator& curPos, const std::vector<double>::const_iterator& endPos);
	virtual void getLocalOptFormat(std::list<W_Point_2>& ptList) const;
	virtual unsigned getNumChildren() const;
	virtual void getTerminalNodeList(std::vector<CSGNode*>& nodeVec);
	
	//! Rotates the points at angles of `theta`, `beta`, and `gamma` about the origin
	void rotatePoints(const Real alpha, const Real beta, const Real gamma);
	//! Translates the points according to the input vector
	void translatePoints(const Vector_3& tp);
	//! Scales points about the origin by an amount of `sf`
	void scalePoints(const Real sf);
	//! Prints debugging information to stdout
	void debugPrint();
	//! Returns the number of points contained here
	int getSize();
	//! Returns priority: A value that sets the overlap order of CSGTerminalNode objects during decoding
	int getPriority() const;
	//! Returns the material properties in a GenericMaterial object
	GenericMaterial getGenericMaterial();
	//! Copies all points into `copyList`
	void copyPointList(std::list<Point_3>& copyList);
	//! Sets the points in this CSGTerminalNode to those in `copyList`
	void setPointList(const std::list<Point_3>& copyList);
	//! Computes the convex hull of the points contained here and puts them in `hullVec`
	void generateHull(Polyhedron_3& P) const;
	//! Computes the convex hull of the points contained here and puts them in `hullVec`
	void generateHull(std::vector<Point_3>& hullVec) const;
	//! If the argument `inpt` point is not contained within the bounding poly, it is projected onto the boundary
	void boundsCheck(Point_3& inpt);
	//! Sets the boundary polygon to `inPoly`
	/*! Points contained in the CSGTerminalNode objects must stay within this bounding polygon */
	static void setBoundPoly(Polyhedron_3& inPoly);
	//! Returns whether or not the argument `testPt` is within the bounding polyhedron
	static bool isPointInBoundary(const Point_3& testPt);
	//! Returns the bounding polyhedron
	static Polyhedron_3 getBoundingPoly();
	//! Returns the number of points that define the bounding polyhedron
	static int getNumBoundingPoints();
	//! Returns the maximum distance between any two points on the bounding poly
	static Real getBoundMag();
	//! Returns a bounding box that contains the bounding polyhedron
	static CGAL::Bbox_3 getBoundBox();

	static int decodecalls;
private:
	// methods
	Point_3 computeHullCenterOMass();
	Point_3 genRandomPoint();
	Point_3 genRandomPointInHull();
	static bool isPointInPoly(const Point_3& tstPt, Polyhedron_3& poly);
	Plane_3 getFacetPlane(Facet_iterator_3& fit) const;
private:
	// data
	static Polyhedron_3 boundingPoly;
	static VerbLevel verbosity;
	std::list<Point_3> geneList;
	GenericMaterial itsGenericMaterial;
	Real priority;
};

inline
void CSGTerminalNode3D::getTerminalNodeList(std::vector<CSGNode*>& nodeVec)
{
	nodeVec.push_back(this);
}

inline
unsigned CSGTerminalNode3D::countNodes() const
{
  return 1;
}

inline
unsigned CSGTerminalNode3D::getNumChildren() const
{
  return 0;
}

inline
int CSGTerminalNode3D::getSize()
{
	return geneList.size();
}

inline
int CSGTerminalNode3D::getPriority() const
{
	return HelperNS::round(priority);
}

inline
GenericMaterial CSGTerminalNode3D::getGenericMaterial()
{
	return itsGenericMaterial;
}

inline
Polyhedron_3 CSGTerminalNode3D::getBoundingPoly()
{
	return boundingPoly;
}

inline
int CSGTerminalNode3D::getNumBoundingPoints()
{
	return boundingPoly.size_of_vertices();
}

inline
Real CSGTerminalNode3D::getBoundMag()
{
	return boundMag;
}
} //namespace
#endif

