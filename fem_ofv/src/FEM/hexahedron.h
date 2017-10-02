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

#ifndef HEXAHEDRON_H
#define HEXAHEDRON_H

#include "cell.h"
#include <cmath>
#include <iostream>
#include <cassert>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "elemedge.h"
#include "element.h"

//! Hexahedron is an abstract class derived class (from Cell) that defines general hexahedra (6-sided polyhedra).
/*! */
class Hexahedron : public Cell
{
public:
	//! Constructor takes a vector of Point3D objects that indicate the Hexahedron's nodes and a GenericMaterial
	/*! As a Hexahedron is not a simplex, node ordering is somewhat tricky.  This assumes that the nodes are not
	 *  in any particular order, though the function to set them to a standard order is not quite right.  In general,
	 *  it is desirable to not have a standard node order enforced outside of the Hexahedron class, so the node order
	 *  processing must be improved.  
	 */
	Hexahedron(CellType inCellType, const std::vector<Point3D*>& inPtVec, const Topologies::GenericMaterial& inMat);
	virtual ~Hexahedron();

	//! @name Data access
	//@{
	virtual Element<Point3D>* getPatch(unsigned short int patchNum) const;
	virtual Point3D* getNode(unsigned short int nodeNum) const;
	virtual unsigned short int getNumNodes() const {return 8;}
	virtual unsigned getNumPatches() const {return 6;}
	virtual bool containsNode(Point3D* inNode) const;
	//@}
	//! @name Element matrix computation functions
	//@{
	virtual Eigen::MatrixXd getElemMat() const = 0;
	virtual Eigen::MatrixXd getLaplacianElemMat(double matVal) const = 0;
	virtual Eigen::MatrixXd getLaplacianElemMat(const Eigen::MatrixXd& matVal) const = 0;
	virtual Eigen::MatrixXd getMassMat() const = 0;
  virtual Eigen::MatrixXd getThermalExpansionMat(double alphaTE) const = 0;
	virtual Eigen::MatrixXd getTFInt() const = 0;
	virtual Eigen::VectorXd getInterpolationVec(double xi, double eta, double zeta) const = 0;
	virtual Eigen::MatrixXd getGradBF() const = 0;
	virtual Eigen::MatrixXd getLaplacianBF1() const = 0;
	virtual Eigen::MatrixXd getLaplacianBF2() const = 0;
	//@}
	//! @name Geometric functions
	//@{
	virtual void justL(double xi1, double xi2, double xi3, Point3D& l1, Point3D& l2, Point3D& l3) const = 0;
	virtual void jacL(double xi1, double xi2, double xi3, double& jac, Point3D& l1, Point3D& l2, Point3D& l3) const = 0;
	virtual Point3D pos(double xi1, double xi2, double xi3) const = 0; 
	virtual double volumeIntegral() const = 0;
	//@}
	//! @name Topological functions
	//@{
	virtual bool isBoundaryCell() const;
	virtual std::unique_ptr<Element<Point3D>> createPatch(unsigned short int patchNum) const;
	virtual void addPatch(Element<Point3D>* thePatch);
	//! Reorders the nodes so that the area coordinate definitions are swapped
	/*! This is to ensure positivity of the Jacobian */
	virtual void switchXi2Xi3();
	virtual void setStandardPatchOrder();
	virtual void addPatchBases(Element<Point3D>* theCPatch, Point3D interpPos[], std::size_t startPoint, int signum);
	virtual void addCellBases(std::size_t& lastAdded);
	virtual void addNodeBases(const Point3D* theNode, std::size_t bfGID);
	//@}
protected:
    // Data Members
    std::vector<Point3D*> nodeVec;
   	std::vector<Element<Point3D>*> patchVec;

	// private functions
	void setInitialNodeOrdering();
	Element<Point3D>* findPatchByNodes(const std::vector<Point3D*>& patchNodeVec) const;
	std::vector<Point3D*> getPatchNodeVec(std::size_t k) const;
	std::vector<std::size_t> getNeighborNodes(std::size_t knode) const;
};

#endif
