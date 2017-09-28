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

#include "mesh3d.h"
#include "lintetra.h"
#include "trilinhex.h"
#include <algorithm>
#include "catch.hpp"

using namespace Topologies;

TEST_CASE( "Testing Mesh3D class with tetraheral elements", "[Mesh3D]" )
{
	// 6 tet mesh
	// From the Delaunay tesselation of a cube
	Point3D p0(0.,0.,0.), p1(1.,0.,0.), p2(0.,1.,0), p3(1.,1.,0.);
	Point3D p4(0.,0.,1.), p5(1.,0.,1.), p6(0.,1.,1), p7(1.,1.,1.);
	std::vector<Point3D*> nodeVec = {&p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7};
	std::vector<double> mats = {1., 2., 3.};
	GenericMaterial gm(mats);
	std::vector<Point3D*> pts1 = {&p4, &p0, &p2, &p1};
	LinTetra tet1(pts1, gm);
	std::vector<Point3D*> pts2 = {&p6, &p4, &p2, &p1};
	LinTetra tet2(pts2, gm);
	std::vector<Point3D*> pts3 = {&p6, &p5, &p4, &p1};
	LinTetra tet3(pts3, gm);
	std::vector<Point3D*> pts4 = {&p6, &p3, &p5, &p1};
	LinTetra tet4(pts4, gm);
	std::vector<Point3D*> pts5 = {&p6, &p2, &p3, &p1};
	LinTetra tet5(pts5, gm);
	std::vector<Point3D*> pts6 = {&p6, &p7, &p5, &p3};
	LinTetra tet6(pts6, gm);
	std::vector<Cell*> elemVec = {&tet1, &tet2, &tet3, &tet4, &tet5, &tet6};
	bool err;
	Mesh3D testMesh3(err, nodeVec, elemVec);
	REQUIRE(!err);
	SECTION("Geometry")
	{
		double tol = 1e-14;
		std::cout << "testMesh3.computeArea(): " << testMesh3.computeArea() << std::endl;

		REQUIRE(fabs(testMesh3.computeArea() - 6.) < tol);
		REQUIRE(fabs(testMesh3.computeVolume() - 1.) < tol);
	}
	SECTION("Topology counts")
	{
		REQUIRE(testMesh3.getNumCells() == 6);
		REQUIRE(testMesh3.getNumPatches() == 18);
		REQUIRE(testMesh3.getNumInternalPatches() == 6);
		REQUIRE(testMesh3.getNumBoundaryPatches() == 12);
		REQUIRE(testMesh3.getNumUnknowns() == 8);
		REQUIRE(testMesh3.getNumElements() == 6);
	}
	SECTION("Basis functions")
	{
		REQUIRE(testMesh3.getGlobalBF(0,0) == 4);
		REQUIRE(testMesh3.getGlobalBF(0,1) == 0);
		REQUIRE(testMesh3.getGlobalBF(0,2) == 2);
		REQUIRE(testMesh3.getGlobalBF(0,3) == 1);
	}
}

TEST_CASE( "Testing Mesh3D class with hexahedral elements", "[Mesh3D]" )
{
	// 4 hex mesh, 1 element in the x-direction, 2 in y & z
	Point3D p000(0.,0.,0.), p010(0.,1.,0.), p020(0.,2.,0.), 
					p100(1.,0.,0.), p110(1.,1.,0.), p120(1.,2.,0.),
					p001(0.,0.,1.), p011(0.,1.,1.), p021(0.,2.,1.),
					p101(1.,0.,1.), p111(1.,1.,1.), p121(1.,2.,1.),
					p002(0.,0.,2.), p012(0.,1.,2.), p022(0.,2.,2.),
					p102(1.,0.,2.), p112(1.,1.,2.), p122(1.,2.,2.);
	std::vector<Point3D*> nodeVec = {&p000, &p010, &p020, &p100, &p110, &p120,
          												 &p001, &p011, &p021, &p101, &p111, &p121, 
																	 &p002, &p012, &p022, &p102, &p112, &p122};
	std::vector<double> mats = {1., 2., 3.};
	GenericMaterial gm(mats);
	std::vector<Point3D*> pts1 = {&p000, &p100, &p110, &p010, &p001, &p101, &p111, &p011};
	TriLinHex hex1(pts1, gm);
	std::vector<Point3D*> pts2 = {&p010, &p110, &p120, &p020, &p011, &p111, &p121, &p021};
	TriLinHex hex2(pts2, gm);
	std::vector<Point3D*> pts3 = {&p001, &p101, &p111, &p011, &p002, &p102, &p112, &p012};
  TriLinHex hex3(pts3, gm);
  std::vector<Point3D*> pts4 = {&p011, &p111, &p121, &p021, &p012, &p112, &p122, &p022};
  TriLinHex hex4(pts4, gm);
	std::vector<Cell*> elemVec = {&hex1, &hex2, &hex3, &hex4};
	bool err;
	Mesh3D testMesh3(err, nodeVec, elemVec);
	REQUIRE(!err);
	SECTION("Geometry")
	{
		double tol = 1e-14;
		REQUIRE(fabs(testMesh3.computeArea() - 16.) < tol);
		REQUIRE(fabs(testMesh3.computeVolume() - 4.) < tol);
	}
	SECTION("Topology counts")
	{
		REQUIRE(testMesh3.getNumCells() == 4);
		REQUIRE(testMesh3.getNumPatches() == 20);
		REQUIRE(testMesh3.getNumInternalPatches() == 4);
		REQUIRE(testMesh3.getNumBoundaryPatches() == 16);
		REQUIRE(testMesh3.getNumUnknowns() == 18);
		REQUIRE(testMesh3.getNumElements() == 4);
	}
	SECTION("Basis functions")
	{
		REQUIRE(testMesh3.getGlobalBF(0,0) == 0);
		REQUIRE(testMesh3.getGlobalBF(0,1) == 3);
		REQUIRE(testMesh3.getGlobalBF(0,2) == 4);
		REQUIRE(testMesh3.getGlobalBF(0,3) == 1);
		REQUIRE(testMesh3.getGlobalBF(0,4) == 6);
		REQUIRE(testMesh3.getGlobalBF(0,5) == 9);
		REQUIRE(testMesh3.getGlobalBF(0,6) == 10);
		REQUIRE(testMesh3.getGlobalBF(0,7) == 7);
	}
}
