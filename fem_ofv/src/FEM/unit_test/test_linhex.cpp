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

#include "trilinhex.h"
#include "point3d.h"
#include "UTIL/genericmaterial.h"
#include <iostream>
#include <fstream>
#include "catch.hpp"

using namespace Topologies;

TEST_CASE( "Testing TriLinHex class", "[TriLinHex]" )
{
	Point3D p0(0.,0.,0.), p1(1.,0.,0.), p2(1.,1.,0.), p3(0.,1.,0.),
					p4(0.,0.,1.), p5(1.,0.,1.), p6(1.,1.,1.), p7(0.,1.,1.);
	std::vector<double> mats(3, 1.);
	GenericMaterial gm(mats);

	SECTION("Starting with correct node ordering")
	{
		std::vector<Point3D*> pts = {&p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7};
		TriLinHex testHex(pts, gm);
		REQUIRE(fabs(testHex.volumeIntegral() - 1.) < 1e-14);
	}
	SECTION("Incorrect node ordering, (negative Jacobian)")
	{
		std::vector<Point3D*> pts = {&p0, &p3, &p2, &p1, &p4, &p7, &p6, &p5};
		TriLinHex testHex(pts, gm);
		REQUIRE(testHex.volumeIntegral() < 0.);
		testHex.switchXi2Xi3();
		REQUIRE(fabs(testHex.volumeIntegral() - 1.) < 1e-14);
	}
	SECTION("Incorrect node ordering, diagonal nodes")
	{
		std::vector<Point3D*> pts = {&p0, &p1, &p6, &p7, &p4, &p3, &p2, &p5};
		TriLinHex testHex(pts, gm);
		REQUIRE(testHex.volumeIntegral() < 0.);
		testHex.switchXi2Xi3();
		REQUIRE(fabs(testHex.volumeIntegral() - 1.) < 1e-14);
  }
	SECTION("Crazy nodes")
	{
		Point3D q0(0.,0.,0.), q1(1.,0.,0.), q2(1.,1.,0.), q3(0.,1.,0.),
						q4(0.,0.,1.), q5(2.,-1.,1.), q6(1.,1.,1.), q7(-1.,2.,1.);
		std::vector<Point3D*> pts = {&q0, &q1, &q2, &q3, &q4, &q5, &q6, &q7};
		TriLinHex testHex(pts, gm);
		REQUIRE(fabs(testHex.volumeIntegral() - 2.) < 1e-14);
	}
// Element matrices
	std::vector<Point3D*> pts = {&p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7};
	TriLinHex testTet(pts, gm);
	SECTION("Stiffness matrix")
	{
		Eigen::MatrixXd elemMat = testTet.getElemMat();
		Eigen::MatrixXd chkMat(24,24);
		std::ifstream matFile("hexelemmat.txt");
		for(unsigned kr = 0; kr < 24; ++kr)
		{
			for(unsigned kc = 0; kc < 24; ++kc)
			{
				double elem;
				matFile >> elem;
				chkMat(kr,kc) = elem;
			}
		}
		matFile.close();
		REQUIRE(elemMat.isApprox(chkMat, 1e-3));
	}
	SECTION("Thermal expansion matrix")
	{
		Eigen::MatrixXd elemMat = testTet.getThermalExpansionMat(1.);
		Eigen::VectorXd dT(8);
		dT << 1.,1.,1.,1.,1.,1.,1.,1.;
		Eigen::VectorXd chkVec(24);
		chkVec.setZero();
		chkVec << 1.25, 1.25, 1.25, -1.25, 1.25, 1.25, -1.25, -1.25, 1.25, 1.25, -1.25, 1.25, 1.25, 1.25, -1.25, -1.25, 1.25, -1.25, -1.25, -1.25, -1.25, 1.25, -1.25, -1.25;
		Eigen::VectorXd res = elemMat*dT;
		std::cout << "res = \n" << res << std::endl;
		REQUIRE(res.isApprox(chkVec,1e-4));
	}
}

