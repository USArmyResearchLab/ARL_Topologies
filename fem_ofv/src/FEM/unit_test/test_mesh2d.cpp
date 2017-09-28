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

#include "mesh2d.h"
#include "lintri.h"
#include "linquad.h"
#include <algorithm>
#include "catch.hpp"

using namespace Topologies;

TEST_CASE( "Testing Mesh2D class with triangular elements", "[Mesh2D]" )
{
	// 2 triangle mesh
	// p3 *------* p2
	//    |tri2 /|
	//    |    / |
	//    |   /  |
	//    |  /   | 
	//    | /    |
	//    |/ tri1|
	// p0 *------* p1
	// Not to scale!  
	Point2D p0(0., 0.), p1(1., 0.), p2(1., 1.), p3(0., 1.);
	std::vector<Point2D*> nodeVec = {&p0, &p1, &p2, &p3};
	std::vector<double> mats = {1., 2., 3.};
	GenericMaterial gm(mats);
	std::vector<Point2D*> pts1 = {nodeVec[0], nodeVec[1], nodeVec[2]};
	LinearTriangle<Point2D> tri1(pts1, gm);
	std::vector<Point2D*> pts2 = {nodeVec[0], nodeVec[2], nodeVec[3]};
	LinearTriangle<Point2D> tri2(pts2, gm);
	std::vector<Element<Point2D>*> elemVec = {&tri1, &tri2};
	bool err;
	Mesh2D testMesh2(err, nodeVec, elemVec);
	REQUIRE(!err);
	SECTION("Geometry")
	{
		double tol = 1e-14;
		REQUIRE(fabs(testMesh2.computeArea() - 1.) < tol);
		REQUIRE(fabs(testMesh2.computeMaterialPropertyTotal(0) - 1.) < tol);
		REQUIRE(fabs(testMesh2.computeMaterialPropertyTotal(1) - 2.) < tol);
		REQUIRE(fabs(testMesh2.computeMaterialPropertyTotal(2) - 3.) < tol);
	}
	SECTION("Topology counts")
	{
		REQUIRE(testMesh2.getNumEdges() == 5);
		REQUIRE(testMesh2.getNumNodes() == 4);
		REQUIRE(testMesh2.getNumInternalEdges() == 1);
		REQUIRE(testMesh2.getNumBoundaryEdges() == 4);
		REQUIRE(testMesh2.getNumUnk() == 4);
		REQUIRE(testMesh2.getNumElements() == 2);
	}
	SECTION("Basis functions")
	{
		std::vector<std::size_t> bfVec(4);
		for(std::size_t k = 0; k < nodeVec.size(); ++k)
			bfVec[k] = testMesh2.getNodeBF(testMesh2.getNode(k));
		std::sort(bfVec.begin(), bfVec.end());
		for(std::size_t k = 0; k < bfVec.size(); ++k)
			REQUIRE(bfVec[k] == k);
	}
}

TEST_CASE( "Testing Mesh2D class with quadrilateral elements", "[Mesh2D]" )
{
	// 4 quad mesh
	//  p02     p12
	//    *------*------* p22
	//    |quad3 |quad4 |
	//    |      |p11   |
	// p01*------*------* p21
	//    |      |      |
	//    |quad1 |quad2 |
	//    *------*------*
	//   p00    p10    p20
	Point2D p00(-1.,-1.), p10(0.,-1.), p20(1.,-1.), 
					p01(-1., 0.), p11(0., 0.), p21(1., 0.),
					p02(-1., 1.), p12(0., 1.), p22(1., 1.);
	std::vector<Point2D*> nodeVec = {&p00, &p10, &p20, &p01, &p11, &p21, &p02, &p12, &p22};
	std::vector<double> mats = {1., 2., 3.};
	GenericMaterial gm(mats);
	std::vector<Point2D*> pts1 = {&p00, &p10, &p11, &p01};
	LinearQuadrilateral<Point2D> quad1(pts1, gm);
	std::vector<Point2D*> pts2 = {&p10, &p20, &p21, &p11};
	LinearQuadrilateral<Point2D> quad2(pts2, gm);
	std::vector<Point2D*> pts3 = {&p01, &p11, &p12, &p02};
	LinearQuadrilateral<Point2D> quad3(pts3, gm);
	std::vector<Point2D*> pts4 = {&p11, &p21, &p22, &p12};
	LinearQuadrilateral<Point2D> quad4(pts4, gm);
	std::vector<Element<Point2D>*> elemVec = {&quad1, &quad2, &quad3, &quad4};
	bool err;
	Mesh2D testMesh2(err, nodeVec, elemVec);
	REQUIRE(!err);
	SECTION("Geometry")
	{
		double tol = 1e-14;
		REQUIRE(fabs(testMesh2.computeArea() - 4.) < tol);
		REQUIRE(fabs(testMesh2.computeMaterialPropertyTotal(0) - 4.) < tol);
		REQUIRE(fabs(testMesh2.computeMaterialPropertyTotal(1) - 2.*4.) < tol);
		REQUIRE(fabs(testMesh2.computeMaterialPropertyTotal(2) - 3.*4.) < tol);
	}
	SECTION("Topology counts")
	{
		REQUIRE(testMesh2.getNumEdges() == 12);
		REQUIRE(testMesh2.getNumNodes() == 9);
		REQUIRE(testMesh2.getNumInternalEdges() == 4);
		REQUIRE(testMesh2.getNumBoundaryEdges() == 8);
		REQUIRE(testMesh2.getNumUnk() == 9);
		REQUIRE(testMesh2.getNumElements() == 4);
	}
	SECTION("Basis functions")
	{
		std::vector<std::size_t> bfVec(9);
		for(std::size_t k = 0; k < nodeVec.size(); ++k)
			bfVec[k] = testMesh2.getNodeBF(testMesh2.getNode(k));
		std::sort(bfVec.begin(), bfVec.end());
		for(std::size_t k = 0; k < bfVec.size(); ++k)
			REQUIRE(bfVec[k] == k);
	}
}
