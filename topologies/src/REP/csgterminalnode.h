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

#ifndef CSGTERMINALNODE_H
#define CSGTERMINALNODE_H

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include "csgnode.h"
#include "genericmaterial.h"
#include "helper.h"

namespace Topologies{
class CSGTerminalNode : public CSGNode
{
public:
	CSGTerminalNode(int curTreeDepth);
	//! Constructor that generates an elliptical set of `numPts` points centered at `center` and with radii `rx` and `ry`
	CSGTerminalNode(int curTreeDepth, unsigned numPts, Point_2_base center, double rx, double ry);
	CSGTerminalNode(const CSGTerminalNode& copy);
	//! A specialized copy constructor that uses a MirrorType to reflect the subtree's nodes in a specified way
	CSGTerminalNode(const CSGTerminalNode& copy, MirrorType inMT);
	//! Constructor that loads a CSGNode from a file
	CSGTerminalNode(int curTreeDepth, std::ifstream& saveFile);
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
	CSGTerminalNode(int curTreeDepth, int numMatParams, int& ptVecPos, int& matPos, int& prPos, int& treeVecPos, 
			  Real* ptArray, Real* matArray, Real* prArray, int* treeArray);
	//! Constructor that sets a subtree up with the specified values, set in vector iterators
	/*! An implementation of the C-style constructor above, uses iterators */
	CSGTerminalNode(int curTreeDepth, int numMatParams, std::vector<double>::const_iterator& ptIt,
                  std::vector<double>::const_iterator& matIt, std::vector<double>::const_iterator& prIt,
                  std::vector<int>::const_iterator& treeIt);
	//! Constructor that copies the points in `ptList`
	CSGTerminalNode(int curTreeDepth, const std::list<W_Point_2>& ptList);
	CSGTerminalNode(CSGTerminalNode && copy);
	CSGTerminalNode& operator=(CSGTerminalNode copy);
	void swap(CSGTerminalNode& arg2);
	virtual ~CSGTerminalNode();

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
	virtual void setVector(std::vector<double>::const_iterator& curPos, const std::vector<double>::const_iterator& endPos);
	virtual void getLocalOptFormat(std::list<W_Point_2>& ptList) const;
	virtual unsigned getNumChildren() const;
	virtual void getTerminalNodeList(std::vector<CSGNode*>& nodeVec);

	//! Prints debugging information to stdout
	void debugPrint();
	//! Returns the number of points contained here
	int getSize();
	//! Returns priority: A value that sets the overlap order of CSGTerminalNode objects during decoding
	int getPriority() const;
	//! Returns the material properties in a GenericMaterial object
	GenericMaterial getGenericMaterial();
	//! Copies all points into `copyList`
	void copyPointList(std::list<W_Point_2>& copyList);
	//! Sets the points in this CSGTerminalNode to those in `copyList`
	void setPointList(const std::list<W_Point_2>& copyList);
	//! Computes the convex hull of the points contained here and puts them in `hullVec`
	void generateHull(std::vector<Point_2>& hullVec) const;
	//! Returns a point at the centroid of the convex hull of the points
	Point_2 computeHullCenterOMass();
	//! Rotates the points an angle of `theta` about the origin
	void rotatePoints(const Real theta);
	//! Translates the points according to the input vector
	void translatePoints(const Vector_2& tp);
	//! Scales points about the origin by an amount of `sf`
	void scalePoints(const Real sf);
	//! Sets the boundary polygon to `inPoly`
	/*! Points contained in the CSGTerminalNode objects must stay within this bounding polygon */
	static void setBoundPoly(std::vector<Point_2> inPoly);
	//! Set whether or not to use the alpha shape of the points as the decoding type, alternative is convex hull
	static void setUseAlpha(bool inUA);
	//! Returns the kth boundary point
	static Point_2 getBoundingPoint(int k);
	//! Returns the number of points in the bounding polygon
	static int getNumBoundingPoints();
	//! Returns the maximum distance between any two points in the bounding polygon
	static Real getBoundMag();

private:
	// methods
	void initMaterials();
	void decode(std::vector<Mesh_Segment_2>& inSegVec) const;
	void decode(std::vector< std::vector<Mesh_Segment_2> >& inSegVec) const;
	bool decodeMaterialFromSegs(Point_2& testPt, const std::vector<Mesh_Segment_2>& segVec, GenericMaterial& theGM) const;
	Nef_polyhedron_2 nefPolyFromSegs(const std::vector<Mesh_Segment_2>& segVec, Nef_polyhedron_2::Boundary boundTreatment) const;
	std::vector<Point_2> getPointVecFromPoly(const std::vector<Mesh_Segment_2>& orderedPoly) const;
	Point_2 wpt2npt(const W_Point_2& inWpt) const;
	Point_2_base npt2wpt(const Point_2& inNpt) const;
	void setAlphaAndWeightsConvex();
	void boundsCheck(Point_2_base& inpt) const;
	void boundsCheck(W_Point_2& inpt) const;
	void boundsCheck(Point_2& inpt) const;
	void initializePoints();
	bool isPtInList(Point_2 inPt, std::vector<Point_2>& ptVec);
	W_Point_2 genRandomPoint() const;
	bool isPointInPoly(Point_2& testPt, std::vector<Point_2>& inPoly) const;
	bool isPointInPoly(Point_2& testPt, std::vector<Mesh_Segment_2>& inPoly) const;
	Real polyArea(std::vector<Mesh_Segment_2>& inPoly) const;
	Real polyPerimeter(std::vector<Mesh_Segment_2>& inPoly) const;
	W_Point_2 constructFixedPoint(Uint kfp) const;
	void snapToBoundary(Point_2& inpt) const;
	Real getCurrentMaxWeight() const;
	std::list<W_Point_2>::iterator getPointIterator(const W_Point_2& findPt);
	void removePointsInList(std::vector<Point_2>& ptVec);

	// data
	static std::vector<Point_2> boundingPoly;
	static Real alphamax, weightmax;
	static bool useAlpha;
	std::list<W_Point_2> geneList;
	GenericMaterial itsGenericMaterial;
	Real priority, alpha;
};

inline
void CSGTerminalNode::getTerminalNodeList(std::vector<CSGNode*>& nodeVec)
{
	nodeVec.push_back(this);
}

inline
unsigned CSGTerminalNode::countNodes() const
{
	return 1;
}

inline
unsigned CSGTerminalNode::getNumChildren() const
{
	return 0;
}

inline
void CSGTerminalNode::getLocalOptFormat(std::list<W_Point_2>& ptList) const
{
	ptList.insert(ptList.end(), geneList.begin(), geneList.end());
}

inline
void CSGTerminalNode::setUseAlpha(bool inUA)
{
	useAlpha = inUA;
}

inline
Point_2_base CSGTerminalNode::npt2wpt(const Point_2& inNpt) const
{
	return Point_2_base(inNpt.x().to_double(), inNpt.y().to_double());
}

inline
Point_2 CSGTerminalNode::wpt2npt(const W_Point_2& inWpt) const
{
	return Point_2(inWpt.x(), inWpt.y());
}

inline
int CSGTerminalNode::getSize()
{
	return geneList.size();
}

inline
int CSGTerminalNode::getPriority() const
{
	return HelperNS::round(priority);
}

inline
GenericMaterial CSGTerminalNode::getGenericMaterial()
{
	return itsGenericMaterial;
}

inline
Point_2 CSGTerminalNode::getBoundingPoint(int k)
{
	return boundingPoly[k];
}

inline
int CSGTerminalNode::getNumBoundingPoints()
{
	return boundingPoly.size();
}

inline
Real CSGTerminalNode::getBoundMag()
{
	return boundMag;
}
} //namespace
#endif

