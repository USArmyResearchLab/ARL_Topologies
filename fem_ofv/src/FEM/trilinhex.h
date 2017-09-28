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

#ifndef TRILINHEX_H
#define TRILINHEX_H

#include "hexahedron.h"
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "UTIL/genericmaterial.h"
#include "point3d.h"

//! LinTetra is a geometrically flat-faced hexahedron, and implements trilinear, nodal basis functions
/*! A trilinear basis function is not strictly linear, but contains higher order terms as follows:
 *  If we have area coordinates \xi_1, \xi_2, and \xi_3, the basis functions are complete to polynomial order 1,
 *  but also contain cross terms \xi_1 \xi_2, \xi_1 \xi_3, \xi_2 \xi_3, and \xi_1 \xi_2 \xi_3
 */
class TriLinHex : public Hexahedron
{
public:
	//! Constructor which takes 8 nodes defining the hexahedron and material properties
	TriLinHex(const std::vector<Point3D*>& inPtVec, const Topologies::GenericMaterial& inMaterial);
	virtual ~TriLinHex(){};

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
	//! Returns the Jacobian at the point (\xi_1, \xi_2, \xi_3)
	double getJacobian(double xi1, double xi2, double xi3) const;
	//@}

private:
	static const unsigned cRho = 0, cLambda = 1, cMu = 2;
	static double xi1Nodes[8], xi2Nodes[8], xi3Nodes[8];
	static double xi1Gauss1[8], xi2Gauss1[8], xi3Gauss1[8], wGauss1[8]; // Gauss rule
	static const unsigned nGauss1 = 8;

	double bf(unsigned k, double xi1, double xi2, double xi3) const {return 0.125*(1 + xi1*xi1Nodes[k])*(1 + xi2*xi2Nodes[k])*(1 + xi3*xi3Nodes[k]);}
	double dbfdxi1(unsigned k, double xi1, double xi2, double xi3) const {return 0.125*xi1Nodes[k]*(1 + xi2*xi2Nodes[k])*(1 + xi3*xi3Nodes[k]);}
	double dbfdxi2(unsigned k, double xi1, double xi2, double xi3) const {return 0.125*xi2Nodes[k]*(1 + xi1*xi1Nodes[k])*(1 + xi3*xi3Nodes[k]);}
	double dbfdxi3(unsigned k, double xi1, double xi2, double xi3) const {return 0.125*xi3Nodes[k]*(1 + xi1*xi1Nodes[k])*(1 + xi2*xi2Nodes[k]);}

	Eigen::MatrixXd dxdxi(double xi1, double xi2, double xi3) const;
	Eigen::MatrixXd dNdxi(double xi1, double xi2, double xi3) const;
	Eigen::MatrixXd elasticityTensor() const;
	Eigen::MatrixXd bmat(double xi1, double xi2, double xi3) const;
	Eigen::MatrixXd dNdx(double xi1, double xi2, double xi3) const;
	Eigen::MatrixXd inverseJac(double xi1, double xi2, double xi3) const;

};
#endif
