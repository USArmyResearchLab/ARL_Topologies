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

#define CATCH_CONFIG_MAIN

#include "lintri.h"
#include "point2d.h"
#include "UTIL/genericmaterial.h"
#include <iostream>
#include "catch.hpp"

using namespace Topologies;

TEST_CASE( "Testing LinearTriangle class, with Point2D", "[LinearTriangle]" )
{
	Point2D p0(0.,0.), p1(1.,0.), p2(0., 1.);
	std::vector<double> mats(3, 1.);
	GenericMaterial gm(mats);

	SECTION("Test correct node ordering")
	{
		std::vector<Point2D*> pts = {&p0, &p1, &p2};
		LinearTriangle<Point2D> testTri(pts, gm);
		REQUIRE(testTri.getArea() == 0.5);
		REQUIRE(testTri.getOrientation() > 0);
	}
	SECTION("Test incorrect node ordering")
	{
		std::vector<Point2D*> pts = {&p0, &p2, &p1};
		LinearTriangle<Point2D> testTri(pts, gm);
		REQUIRE(testTri.getArea() == 0.5);
		REQUIRE(testTri.getOrientation() < 0);
		testTri.switchXiEta();
		REQUIRE(testTri.getOrientation() > 0);
	}
	// Element matrices
	std::vector<Point2D*> pts = {&p0, &p1, &p2};
	LinearTriangle<Point2D> testTri(pts, gm);
	SECTION("Stiffness matrix")
	{
		Eigen::MatrixXd elemMat = testTri.getElemMat(eptPlaneStrain);
		// Element matrix taken from a Matlab FEM code
		Eigen::MatrixXd chckMat(6,6);
		chckMat <<   2.0, -1.5, -0.5,  1.0, -0.5, -0.5,
								-1.5,  1.5,  0.0, -0.5,  0.0,  0.5,
								-0.5,  0.0,  0.5, -0.5,  0.5,  0.0,
								 1.0, -0.5, -0.5,  2.0, -0.5, -1.5,
								 -0.5,  0.0,  0.5, -0.5,  0.5,  0.0,
								 -0.5,  0.5,  0.0, -1.5,  0.0,  1.5;
		REQUIRE(elemMat.isApprox(chckMat, 1e-14));
	}
	SECTION("Thermal expansion matrix")
	{
		Eigen::MatrixXd elemMat = testTri.getThermalExpansionMat(eptPlaneStrain, 1.);
		Eigen::VectorXd dT(3);
		dT << 1.,1.,1.;
		Eigen::VectorXd chkVec(6);
		chkVec << 2., -2., 0., 2., 0., -2.;
		Eigen::VectorXd res = elemMat*dT;
		REQUIRE(res.isApprox(chkVec,1e-4));
	}
}
