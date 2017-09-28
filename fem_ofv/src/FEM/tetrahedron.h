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

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

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

//! Tetrahedron is an abstract derived class (from Cell), which defines  geometry and numerics for Tetrahedral elements
/*! Currently, Tetrahedron has one derived class, LinTetra, which implements linear, nodal basis functions.  In the future,
 * this class could be extended to include higher order elements.
 */
class Tetrahedron : public Cell
{
public:
	//! Constructor that takes a vector of four points, defining the tetrahedron and a GenericMaterial with material properties for the Tetrahedron
	/*! In terms of normalized area coordinates, the point inPtVec[0] is the \xi_4 = 1 point, inPtVec[1] is the \xi_1 = 1 point, 
	 *  inPtVec[2] is the \xi_2 = 1 point, and inPtVec[3] is the \xi_3 = 1 point. 
	 */
	Tetrahedron(const std::vector<Point3D*>& inPtVec, const Topologies::GenericMaterial& inMat);
	virtual ~Tetrahedron();

	//! @name Data access
	//@{
	virtual unsigned short int getNumNodes() const {return 4;}
	virtual unsigned getNumPatches() const {return 4;}
	virtual Element<Point3D>* getPatch(unsigned short int patchNum) const;
	virtual Point3D* getNode(unsigned short int nodeNum) const;
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
	//! Switches the node associated with the point xi2 = 1 with the node at xi3 = 1
	virtual void switchXi2Xi3();
	//! Reorder patches so they are listed in the order of the node opposite to them. Once this is done, patchVec[0] is the patch opposite nodeVec[0], etc.
	virtual void setStandardPatchOrder();
	virtual void addPatchBases(Element<Point3D>* theCPatch, Point3D interpPos[], std::size_t startPoint, int signum);
	virtual void addCellBases(std::size_t& lastAdded);
	virtual void addNodeBases(const Point3D* theNode, std::size_t bfGID);
	//@}

protected:
	// Data Members
	std::vector<Point3D*> nodeVec;
	std::vector<Element<Point3D>*> patchVec;
};

#endif
//
