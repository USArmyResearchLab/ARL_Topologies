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

#include "boundarycondition.h"
#include "REP/tomesh.h"
#include "IO/exotxtmeshloader.h"
#include "IO/gmshtxtloader.h"
#include "catch.hpp"
#include <iostream>

using namespace Topologies;

TEST_CASE("Testing BoundaryCondition on a 2 element tri mesh", "[BoundaryCondition]")
{
	// Mesh
	Point_2_base p0(0., 0.), p1(1., 0.), p2(1., 1.), p3(0., 1.);
	std::vector<Point_2_base> nodeVec = {p0, p1, p2, p3};
	std::vector<std::vector<std::size_t>> connVec(2);
	std::vector<std::size_t> cv1 = {0, 1, 2};
	connVec[0] = cv1;
	std::vector<std::size_t> cv2 = {0, 2, 3};
	connVec[1] = cv2;
	TOMesh2D testMesh(nodeVec, connVec);
	SECTION("Vertical line at 1")
	{
		// Geometry of boundary condition
		std::unique_ptr<GeometricEntity> upGE(new InfiniteLine(1., false));
		BoundaryCondition testBC(bctVLine, true, true, std::move(upGE));
		REQUIRE(testBC.checkValidity(&testMesh));
		std::vector<std::size_t> bcPts = testBC.applyBC(&testMesh);
		REQUIRE(bcPts.size() == 2);
		REQUIRE(bcPts[0] == 1);
		REQUIRE(bcPts[1] == 2);
	}
	SECTION("Vertical line at 2")
	{
		// Geometry of boundary condition
		std::unique_ptr<GeometricEntity> upGE(new InfiniteLine(2., false));
		BoundaryCondition testBC(bctVLine, true, true, std::move(upGE));
		REQUIRE(!testBC.checkValidity(&testMesh));
		std::vector<std::size_t> bcPts = testBC.applyBC(&testMesh);
		REQUIRE(bcPts.empty());
	}
}

TEST_CASE("Testing 2D GMSH mesh", "[BoundaryCondition]")
{
	std::unique_ptr<TOMesh2D> upTestMesh = InputLoader::GMSHTxtMeshLoader::loadGMSHTxtFile2D("testtri.msh");
	BoundaryCondition testBC(bctUnknown, true, true, true, 8, mffGMSH, "testtri.msh", 2);
	REQUIRE(testBC.checkValidity(upTestMesh.get()));
	std::vector<std::size_t> bcPts = testBC.applyBC(upTestMesh.get());
	REQUIRE(bcPts.size() == 11);
	REQUIRE(bcPts[0] == 0);
	REQUIRE(bcPts[1] == 1);
	REQUIRE(bcPts[2] == 31);
	REQUIRE(bcPts[3] == 32);
	REQUIRE(bcPts[4] == 33);
	REQUIRE(bcPts[5] == 34);
	REQUIRE(bcPts[6] == 35);
	REQUIRE(bcPts[7] == 36);
	REQUIRE(bcPts[8] == 37);
	REQUIRE(bcPts[9] == 38);
	REQUIRE(bcPts[10] == 39);
}

TEST_CASE("Testing 2D Exodus mesh", "[BoundaryCondition]")
{
	std::unique_ptr<TOMesh2D> upTestMesh = InputLoader::ExoTxtMeshLoader::loadExoTxtFileTri("bctesttrimesh.txt");
	BoundaryCondition testBC(bctUnknown, true, true, true, 1, mffExodus, "bctesttrimesh.txt", 2);
	REQUIRE(testBC.checkValidity(upTestMesh.get()));
	std::vector<std::size_t> bcPts = testBC.applyBC(upTestMesh.get());
	REQUIRE(bcPts.size() == 2);
	REQUIRE(bcPts[0] == 0);
	REQUIRE(bcPts[1] == 1);
}

TEST_CASE("Testing BoundaryCondition on a 1 element quad mesh", "[BoundaryCondition]")
{
	// Mesh
	Point_2_base p0(0., 0.), p1(1., 0.), p2(1., 1.), p3(0., 1.);
	std::vector<Point_2_base> nodeVec = {p0, p1, p2, p3};
	std::vector<std::vector<std::size_t>> connVec(1);
	std::vector<std::size_t> cv1 = {0, 1, 2, 3};
	connVec[0] = cv1;
	TOMesh2D testMesh(nodeVec, connVec);
	SECTION("Vertical line at 1")
	{
		// Geometry of boundary condition
		std::unique_ptr<GeometricEntity> upGE(new InfiniteLine(1., false));
		BoundaryCondition testBC(bctVLine, true, true, std::move(upGE));
		REQUIRE(testBC.checkValidity(&testMesh));
		std::vector<std::size_t> bcPts = testBC.applyBC(&testMesh);
		REQUIRE(bcPts.size() == 2);
		REQUIRE(bcPts[0] == 1);
		REQUIRE(bcPts[1] == 2);
	}
	SECTION("Vertical line at 2")
	{
		// Geometry of boundary condition
		std::unique_ptr<GeometricEntity> upGE(new InfiniteLine(2., false));
		BoundaryCondition testBC(bctVLine, true, true, std::move(upGE));
		REQUIRE(!testBC.checkValidity(&testMesh));
		std::vector<std::size_t> bcPts = testBC.applyBC(&testMesh);
		REQUIRE(bcPts.empty());
	}
}

TEST_CASE("Testing BoundaryCondition on a 1 element tet mesh", "[BoundaryCondition]")
{
	std::cout << "Testing 3d tet mesh" << std::endl;
  // Mesh
  Point_3_base p0(0.,0.,0.), p1(1.,0.,0.), p2(0.,1.,0.), p3(0.,0.,1.);
  std::vector<Point_3_base> nodeVec = {p0, p1, p2, p3};
  std::vector<std::vector<std::size_t>> connVec(1);
  std::vector<std::size_t> cv1 = {0, 1, 2, 3};
  connVec[0] = cv1;
  TOMesh3D testMesh(nodeVec, connVec);
  SECTION("Infinite xy plane at z=0")
  {
    // Geometry of boundary condition
		std::unique_ptr<GeometricEntity> upGE(new InfinitePlane(0., GeometryEntityFactory::poXY));
//    std::unique_ptr<GeometricEntity> upGE(new InfinitePlane(p0, p3));
    BoundaryCondition testBC(bctInfinitePlane, true, true, true, std::move(upGE));
    REQUIRE(testBC.checkValidity(&testMesh));
    std::vector<std::size_t> bcPts = testBC.applyBC(&testMesh);
    REQUIRE(bcPts.size() == 3);
    REQUIRE(bcPts[0] == 0);
    REQUIRE(bcPts[1] == 1);
		REQUIRE(bcPts[2] == 2);
  }
  SECTION("Infinite xy plane at z=2")
  {
    // Geometry of boundary condition
		Point_3_base q(0.,0.,2.);
    std::unique_ptr<GeometricEntity> upGE(new InfinitePlane(q, p3));
    BoundaryCondition testBC(bctInfinitePlane, true, true, true, std::move(upGE));
    REQUIRE(!testBC.checkValidity(&testMesh));
    std::vector<std::size_t> bcPts = testBC.applyBC(&testMesh);
    REQUIRE(bcPts.empty());
  }
}
