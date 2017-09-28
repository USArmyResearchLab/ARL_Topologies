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

#include "lintetra.h"
#include "point3d.h"
#include "UTIL/genericmaterial.h"
#include <iostream>
#include "catch.hpp"

using namespace Topologies;

TEST_CASE( "Testing LinTetra class", "[LinTetra]" )
{
	Point3D p0(0.,0.,0.), p1(1.,0.,0.), p2(0.,1.,0.), p3(0.,0.,1.);
	std::vector<double> mats(3, 1.);
	GenericMaterial gm(mats);

	SECTION("Starting with correct node ordering")
	{
		std::vector<Point3D*> pts = {&p0, &p1, &p2, &p3};
		LinTetra testTetra(pts, gm);
		REQUIRE(fabs(testTetra.volumeIntegral() - 1./6.) < 1e-14);
	}
	SECTION("Incorrect node ordering, (negative Jacobian)")
	{
		std::vector<Point3D*> pts = {&p0, &p3, &p2, &p1};
		LinTetra testTetra(pts, gm);
		REQUIRE(testTetra.volumeIntegral() < 0.);
		testTetra.switchXi2Xi3();
		REQUIRE(fabs(testTetra.volumeIntegral() - 1./6.) < 1e-14);
	}
// Element matrices
  std::vector<Point3D*> pts = {&p0, &p1, &p2, &p3};
  LinTetra testTet(pts, gm);
  SECTION("Stiffness matrix")
  {
    Eigen::MatrixXd elemMat = testTet.getElemMat();
		Eigen::MatrixXd chkMat(12,12);
		chkMat << 0.8333, 0.3333, 0.3333, -0.5000, -0.1667, -0.1667, -0.1667, -0.1667,      0, -0.1667,  0, -0.1667,
    0.3333,   0.8333,   0.3333,  -0.1667,  -0.1667,        0,  -0.1667,  -0.5000,  -0.1667,        0,  -0.1667,  -0.1667,
    0.3333,   0.3333,   0.8333,  -0.1667,        0,  -0.1667,        0,   -0.1667,   -0.1667,   -0.1667,   -0.1667,   -0.5000,
   -0.5000,  -0.1667,   -0.1667,    0.5000,         0,         0,         0,    0.1667,         0,         0,         0,   0.1667,
   -0.1667,  -0.1667,         0,         0,    0.1667,         0,    0.1667,         0,         0,         0,         0,        0,
   -0.1667,        0,   -0.1667,         0,         0,    0.1667,         0,         0,         0,    0.1667,         0,        0,
   -0.1667,  -0.1667,         0,         0,    0.1667,         0,    0.1667,         0,         0,         0,         0,        0,
   -0.1667,  -0.5000,   -0.1667,    0.1667,         0,         0,         0,    0.5000,         0,         0,         0,   0.1667,
         0,  -0.1667,   -0.1667,         0,         0,         0,         0,         0,    0.1667,         0,    0.1667,        0,
   -0.1667,        0,   -0.1667,         0,         0,    0.1667,         0,         0,         0,    0.1667,         0,        0,
         0,  -0.1667,   -0.1667,         0,         0,         0,         0,         0,    0.1667,         0,    0.1667,        0,
   -0.1667,  -0.1667,   -0.5000,    0.1667,         0,         0,         0,    0.1667,         0,         0,         0,   0.5000;
		REQUIRE(elemMat.isApprox(chkMat, 1e-3));
	}
	SECTION("Thermal expansion matrix")
  {
    Eigen::MatrixXd elemMat = testTet.getThermalExpansionMat(1.);
    Eigen::VectorXd dT(4);
    dT << 1.,1.,1.,1.;
    Eigen::VectorXd chkVec(12);
    chkVec.setZero();
    chkVec << 5./6., 5./6., 5./6., -5./6., 0., 0., 0., -5./6., 0., 0., 0., -5./6.;
    Eigen::VectorXd res = elemMat*dT;
    std::cout << "res = \n" << res << std::endl;
    REQUIRE(res.isApprox(chkVec,1e-4));
  }
}
