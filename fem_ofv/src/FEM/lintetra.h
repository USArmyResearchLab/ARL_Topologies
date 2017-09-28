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

#ifndef LINTETRA_H
#define LINTETRA_H

#include "tetrahedron.h"
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "genericmaterial.h"
#include "point3d.h"

//! LinTetra is a geometrically flat-faced tetrahedron, and implements linear, nodal basis functions
class LinTetra : public Tetrahedron
{
public:
    //! Constructor that takes a vector of four points, defining the tetrahedron and a GenericMaterial with material properties for the Tetrahedron
	/*! In terms of normalized area coordinates, the point inPtVec[0] is the \xi_4 = 1 point, inPtVec[1] is the \xi_1 = 1 point, 
	 *  inPtVec[2] is the \xi_2 = 1 point, and inPtVec[3] is the \xi_3 = 1 point. 
	 */
	LinTetra(const std::vector<Point3D*>& inPtVec, const Topologies::GenericMaterial& inMaterial);
	virtual ~LinTetra(){};

    //! @name Element matrix computation functions
	//@{
	virtual Eigen::MatrixXd getElemMat() const;
	virtual Eigen::MatrixXd getLaplacianElemMat(double matVal) const;
	virtual Eigen::MatrixXd getLaplacianElemMat(const Eigen::MatrixXd& matVal) const;
	virtual Eigen::MatrixXd getMassMat() const;
	virtual Eigen::MatrixXd getThermalExpansionMat(double alphaTE) const;
	virtual Eigen::MatrixXd getTFInt() const;
	virtual Eigen::VectorXd getInterpolationVec(double xi, double eta, double zeta) const;
	virtual Eigen::MatrixXd getGradBF() const;
	virtual Eigen::MatrixXd getLaplacianBF1() const;
	virtual Eigen::MatrixXd getLaplacianBF2() const;
	//@}
	//! @name Geometric functions
	//@{
	virtual void justL(double xi1, double xi2, double xi3, Point3D& l1, Point3D& l2, Point3D& l3) const;
	virtual void jacL(double xi1, double xi2, double xi3, double& jac, Point3D& l1, Point3D& l2, Point3D& l3) const;
	virtual void switchXi2Xi3();
	virtual Point3D pos(double xi1, double xi2, double xi3) const;
	virtual double volumeIntegral() const;
	//! Returns the element Jacobian.  Coordinates are not needed as the Jacobian is constant over the entire tetrahedron.
	double getJacobian() const;
	//@}
private:
	Eigen::MatrixXd dxdxi() const;
	Eigen::MatrixXd dNdxi() const;
	Eigen::MatrixXd elasticityTensor() const;
	Eigen::MatrixXd bmat() const;
	Eigen::MatrixXd dNdx() const;
	Eigen::MatrixXd inverseJac() const;
	// Data Members
	Point3D r, rXi1, rXi2, rXi3;
};
#endif
