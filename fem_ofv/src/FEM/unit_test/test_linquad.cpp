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

#include "linquad.h"
#include "point2d.h"
#include "UTIL/genericmaterial.h"
#include <iostream>
#include "catch.hpp"

using namespace Topologies;

TEST_CASE( "Testing LinearQuadrilateral class, with Point2D", "[LinearQuadrilateral]" )
{
	Point2D p0(0.,0.), p1(1.,0.), p2(1., 1.), p3(0., 1.);
	std::vector<double> mats(3, 1.);
	GenericMaterial gm(mats);

	SECTION("Starting with correct node ordering")
	{
		std::vector<Point2D*> pts = {&p0, &p1, &p2, &p3};
		LinearQuadrilateral<Point2D> testQuad(pts, gm);
		REQUIRE(testQuad.getArea() == 1.);
		REQUIRE(testQuad.getOrientation() > 0);
	}
	SECTION("Opposite node ordering, (negative Jacobian)")
	{
		std::vector<Point2D*> pts = {&p0, &p3, &p2, &p1};
		LinearQuadrilateral<Point2D> testQuad(pts, gm);
		REQUIRE(fabs(testQuad.getArea()) == 1.);
		REQUIRE(testQuad.getOrientation() < 0);
		testQuad.switchXiEta();
		REQUIRE(testQuad.getOrientation() > 0);
	}
	SECTION("Incorrect node ordering (intersecting segments")
  {
		std::vector<Point2D*> pts = {&p0, &p2, &p3, &p1};
		LinearQuadrilateral<Point2D> testQuad(pts, gm);
		REQUIRE(fabs(testQuad.getArea() == 1.));
	}
	// Element matrices
	std::vector<Point2D*> pts = {&p0, &p1, &p2, &p3};
	LinearQuadrilateral<Point2D> testQuad(pts, gm);
	SECTION("Stiffness matrix")
	{
		Eigen::MatrixXd elemMat = testQuad.getElemMat(eptPlaneStrain);
		Eigen::MatrixXd chkMat(8,8);
		chkMat << 1.3333,-0.8333,-0.6667, 0.1667, 0.5000, 0.0000,-0.5000, 0.000,
					   -0.8333, 1.3333, 0.1667,-0.6667, 0.0000,-0.5000,-0.0000, 0.5000,
						 -0.6667, 0.1667, 1.3333,-0.8333,-0.5000, 0.0000, 0.5000, 0.0000,
							0.1667,-0.6667,-0.8333, 1.3333, 0.0000, 0.5000, 0.0000,-0.5000,
							0.5000, 0.0000,-0.5000, 0.0000, 1.3333, 0.1667,-0.6667,-0.8333,
							0.0000,-0.5000,-0.0000, 0.5000, 0.1667, 1.3333,-0.8333,-0.6667,
						 -0.5000, 0.0000, 0.5000, 0.0000,-0.6667,-0.8333, 1.3333, 0.1667,
							0.0000, 0.5000, 0.0000,-0.5000,-0.8333,-0.6667, 0.1667, 1.3333;
		REQUIRE(elemMat.isApprox(chkMat, 1e-4));
	}
	SECTION("Thermal expansion matrix")
	{
		Eigen::MatrixXd elemMat = testQuad.getThermalExpansionMat(eptPlaneStrain, 1.);
		Eigen::VectorXd dT(4);
		dT << 1,1,1,1;
		Eigen::VectorXd chkVec(8);
		chkVec << 2., -2., -2., 2., 2., 2., -2., -2.;
		Eigen::VectorXd res = elemMat*dT;
		REQUIRE(res.isApprox(chkVec,1e-4));
	}
}
