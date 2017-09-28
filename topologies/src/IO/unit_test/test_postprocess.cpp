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

// Questions?
// Contact: Raymond Wildman, raymond.a.wildman.civ@mail.mil

#define CATCH_CONFIG_MAIN

#include "postprocess.h"
#include "tomesh.h"
#include "cartesianmesher.h"
#include "catch.hpp"

using namespace Topologies;

TEST_CASE("Testing 2d simplification","[PostProcess]")
{
	using namespace PostProcess;
	Point_2_base p0(-1.,-1.), p1(0.,0.), p2(1.,1.), p3(1.,0.);
	std::vector<Mesh_Segment_2> inVec, outVec;
	SECTION("2 segment simplification")
	{
		inVec = {Mesh_Segment_2(p0, p1), Mesh_Segment_2(p1, p2)};
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 1);
		REQUIRE(outVec[0].source().x() == Approx(p0.x()));
		REQUIRE(outVec[0].source().y() == Approx(p0.y()));
		REQUIRE(outVec[0].target().x() == Approx(p2.x()));
		REQUIRE(outVec[0].target().y() == Approx(p2.y()));
		// Add non-parallel segment
		inVec.push_back(Mesh_Segment_2(p2, p3));
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 2);
		REQUIRE(outVec[0].source().x() == Approx(p0.x()));
		REQUIRE(outVec[0].source().y() == Approx(p0.y()));
		REQUIRE(outVec[0].target().x() == Approx(p2.x()));
		REQUIRE(outVec[0].target().y() == Approx(p2.y()));
		REQUIRE(outVec[1].target().x() == Approx(p3.x()));
		REQUIRE(outVec[1].target().y() == Approx(p3.y()));
		// Change order
		inVec = {Mesh_Segment_2(p3, p2), Mesh_Segment_2(p2, p1), Mesh_Segment_2(p1, p0)};
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 2);
		REQUIRE(outVec[0].source().x() == Approx(p3.x()));
		REQUIRE(outVec[0].source().y() == Approx(p3.y()));
		REQUIRE(outVec[1].target().x() == Approx(p0.x()));
		REQUIRE(outVec[1].target().y() == Approx(p0.y()));
		// Cyclic, simplified edge fist
		inVec = {Mesh_Segment_2(p0, p1), Mesh_Segment_2(p1, p2), Mesh_Segment_2(p2, p3), Mesh_Segment_2(p3, p0)};
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 3);
		// Cyclic, simplified edge split
		inVec = {Mesh_Segment_2(p1, p2), Mesh_Segment_2(p2, p3), Mesh_Segment_2(p3, p0), Mesh_Segment_2(p0, p1)};
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 3);
	}
	SECTION("No simplification")
	{
		inVec = {Mesh_Segment_2(p0, p1), Mesh_Segment_2(p1, p3)};
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 2);
		// Cyclic
		inVec = {Mesh_Segment_2(p0, p1), Mesh_Segment_2(p1, p3), Mesh_Segment_2(p3, p0)};
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 3);
	}
	SECTION("Check tolerance")
	{
		double tol = 1e-6;
		p3 = Point_2_base(1., 1. + tol);
		inVec = {Mesh_Segment_2(p0, p1), Mesh_Segment_2(p1, p3)};
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 1);
		tol = 1e-5;
		p3 = Point_2_base(1., 1. + tol);
		inVec = {Mesh_Segment_2(p0, p1), Mesh_Segment_2(p1, p3)};
		simplify2d(inVec, outVec);
		REQUIRE(outVec.size() == 2);
	}
}

TEST_CASE("Testing 2d subdivision","[PostProcess]")
{
	// Check that boundary segs are not subdivided
	using namespace PostProcess;
	Point_2_base p0(-1.,-1.), p1(1.,-1.), p2(1.,1.), p3(-1.,1.);
	Point_2_base q0(0.,-1.), q1(0., -0.5), q2(0.,0.), q3(-0.5, 0.), q4(-1.,0.);
	// Check no change in coarse (all on boundary)
	std::vector<Mesh_Segment_2> coarse = {Mesh_Segment_2(p0, q0), Mesh_Segment_2(q0, q2), 
		Mesh_Segment_2(q2, q4), Mesh_Segment_2(q4, p0)};
	std::vector<Mesh_Segment_2> boundary = {Mesh_Segment_2(p0, p1), Mesh_Segment_2(p1, p2),
		Mesh_Segment_2(p2, p3), Mesh_Segment_2(p3, p0)};
	std::vector<Mesh_Segment_2> fine;
	subdivision2d(coarse, fine, boundary);
	REQUIRE(fine.size() == 4);
	// Check with subdivision
	coarse = {Mesh_Segment_2(p0, q0), Mesh_Segment_2(q0, q1), Mesh_Segment_2(q1, q2),
		Mesh_Segment_2(q2, q3), Mesh_Segment_2(q3, q4), Mesh_Segment_2(q4, p0)};
	subdivision2d(coarse, fine, boundary);
	REQUIRE(fine.size() == 8);
	REQUIRE(coarse[0].source().x() == Approx(fine[0].source().x()));
	REQUIRE(coarse[0].source().y() == Approx(fine[0].source().y()));
	REQUIRE(coarse[0].target().x() == Approx(fine[0].target().x()));
	REQUIRE(coarse[0].target().y() == Approx(fine[0].target().y()));
	REQUIRE(coarse[1].source().x() == Approx(fine[1].source().x()));
	REQUIRE(coarse[1].source().y() == Approx(fine[1].source().y()));
	REQUIRE(coarse[1].target().x() == Approx(fine[1].target().x()));
	REQUIRE(coarse[1].target().y() == Approx(fine[1].target().y()));
	REQUIRE(coarse[2].source().x() == Approx(fine[2].source().x()));
	REQUIRE(coarse[2].source().y() == Approx(fine[2].source().y()));
	REQUIRE(coarse[2].target().x() == Approx(fine[2].target().x()));
	REQUIRE_FALSE(coarse[2].target().y() == Approx(fine[2].target().y()));
	REQUIRE_FALSE(coarse[3].source().x() == Approx(fine[5].source().x()));
	REQUIRE(coarse[3].source().y() == Approx(fine[5].source().y()));
	REQUIRE(coarse[3].target().x() == Approx(fine[5].target().x()));
	REQUIRE(coarse[3].target().y() == Approx(fine[5].target().y()));
	REQUIRE(coarse[4].source().x() == Approx(fine[6].source().x()));
	REQUIRE(coarse[4].source().y() == Approx(fine[6].source().y()));
	REQUIRE(coarse[4].target().x() == Approx(fine[6].target().x()));
	REQUIRE(coarse[4].target().y() == Approx(fine[6].target().y()));
	REQUIRE(coarse[5].source().x() == Approx(fine[7].source().x()));
	REQUIRE(coarse[5].source().y() == Approx(fine[7].source().y()));
	REQUIRE(coarse[5].target().x() == Approx(fine[7].target().x()));
	REQUIRE(coarse[5].target().y() == Approx(fine[7].target().y()));
}

TEST_CASE("Testing 2d iso-surface generation","[PostProcess]")
{
	using namespace PostProcess;
	Point_2_base p0(-1.,-1.), p1(1.,-1.), p2(1.,1.), p3(-1.,1.);
	std::vector<Point_2_base> ptVec = {p0, p1, p2, p3};
	std::vector<double> valVec(4, 1.);
	// Test all below threshold
	std::vector<Mesh_Segment_2> res = meshSegsIsoSurf2d(valVec, ptVec, 2.);
	REQUIRE(res.empty());
	// Test all above
	res = meshSegsIsoSurf2d(valVec, ptVec, 0.);
	REQUIRE(res.empty());
	// Test line at x=0
	valVec = {0., 2., 2., 0.};
	res = meshSegsIsoSurf2d(valVec, ptVec, 1.);
	REQUIRE_FALSE(res.empty());
	for(auto it = res.begin(); it != res.end(); ++it)
	{
		REQUIRE(it->source().x() == Approx(0.));
		REQUIRE(it->target().x() == Approx(0.));
	}
	// Test line at y=0
	valVec = {0., 0., 2., 2.};
	res = meshSegsIsoSurf2d(valVec, ptVec, 1.);
	REQUIRE_FALSE(res.empty());
	for(auto it = res.begin(); it != res.end(); ++it)
	{
		REQUIRE(it->source().y() == Approx(0.));
		REQUIRE(it->target().y() == Approx(0.));
	}
}

TEST_CASE("Testing 3d interpolation","[PostProcess]")
{
	using namespace PostProcess;
	// Test 1 element mesh
	std::vector<double> elemVals(1, 1.);
	std::unique_ptr<TOMesh3D> testMesh = CartesianMesher::generateMesh(metHex, 1, 1, 1, 1., 1., 1., elemVals);
	std::vector<double> nodalVals(testMesh->getNumNodes(), 1.);
	TetMeshInterp testTMI(nodalVals, testMesh.get(), 0.5);
	Point_3_base p0(0.0001,0.0001,0.0001), p1(0.5,0.5,0.5), p2(0.9999, 0.9999, 0.9999), p3(1.001, 1.001, 1.001);
	REQUIRE(testTMI(p0) == Approx(-0.5));
	REQUIRE(testTMI(p1) == Approx(-0.5));
	REQUIRE(testTMI(p2) == Approx(-0.5));
	REQUIRE(testTMI(p3) == Approx(1.0));
	// 1x3x1
	elemVals = std::vector<double>(3, 1.);
	testMesh = CartesianMesher::generateMesh(metHex, 1, 3, 1, 1., 1., 1., elemVals);
	nodalVals = std::vector<double>(testMesh->getNumNodes(), 1.);
	testTMI = TetMeshInterp(nodalVals, testMesh.get(), 0.5);
	REQUIRE(testTMI(p0) == Approx(-0.5));
	REQUIRE(testTMI(p1) == Approx(-0.5));
	REQUIRE(testTMI(p2) == Approx(-0.5));
	REQUIRE(testTMI(p3) == Approx(1.0));
	// 2x2x1
	elemVals = std::vector<double>(4, 1.);
	testMesh = CartesianMesher::generateMesh(metHex, 2, 2, 1, 1., 1., 1., elemVals);
	nodalVals = std::vector<double>(testMesh->getNumNodes(), 1.);
	testTMI = TetMeshInterp(nodalVals, testMesh.get(), 0.5);
	REQUIRE(testTMI(p0) == Approx(-0.5));
	REQUIRE(testTMI(p1) == Approx(-0.5));
	REQUIRE(testTMI(p2) == Approx(-0.5));
	REQUIRE(testTMI(p3) == Approx(1.0));
	// 2x2x2
	elemVals = std::vector<double>(8, 1.);
	testMesh = CartesianMesher::generateMesh(metHex, 2, 2, 2, 1., 1., 1., elemVals);
	nodalVals = std::vector<double>(testMesh->getNumNodes(), 1.);
	testTMI = TetMeshInterp(nodalVals, testMesh.get(), 0.5);
	REQUIRE(testTMI(p0) == Approx(-0.5));
	REQUIRE(testTMI(p1) == Approx(-0.5));
	REQUIRE(testTMI(p2) == Approx(-0.5));
	REQUIRE(testTMI(p3) == Approx(1.0));
	// 1x1x2
	elemVals = std::vector<double>(2, 1.);
	testMesh = CartesianMesher::generateMesh(metHex, 1, 1, 2, 1., 1., 2., elemVals);
	nodalVals = {0., 0., 0., 0., 0.5, 0.5, 0.5, 0.5, 1., 1., 1., 1.};
	testTMI = TetMeshInterp(nodalVals, testMesh.get(), 0.5);
	Point_3_base q0(0.5, 0.5, 0.1), q1(0.5, 0.5, 1.), q2(0.5, 0.5, 1.9);
	REQUIRE(testTMI(q0) == Approx( 0.45));
	REQUIRE(testTMI(q1) == Approx(0.));
	REQUIRE(testTMI(q2) == Approx(-0.45));
	// Repeat with tetrahedral elements
	testMesh = CartesianMesher::generateMesh(metTet, 1, 1, 2, 1., 1., 2., elemVals);
	nodalVals = {0., 0., 0., 0., 0.5, 0.5, 0.5, 0.5, 1., 1., 1., 1.};
	testTMI = TetMeshInterp(nodalVals, testMesh.get(), 0.5);
	REQUIRE(testTMI(q0) == Approx( 0.45));
	REQUIRE(testTMI(q1) == Approx(0.));
	REQUIRE(testTMI(q2) == Approx(-0.45));
	// TODO: Add a test with a concave mesh
}

TEST_CASE("Testing 3d voxel interpolation","[PostProcess]")
{
  using namespace PostProcess;
  // Test 1 element mesh
  std::vector<double> elemVals(1, 1.);
  std::vector<double> nodalVals(8, 1.);
  VoxelInterp testTMI(nodalVals, 1, 1, 1, 1., 1., 1., 0.5);
  Point_3_base p0(0.0001,0.0001,0.0001), p1(0.5,0.5,0.5), p2(0.9999, 0.9999, 0.9999), p3(1.001, 1.001, 1.001);
  REQUIRE(testTMI(p0) == Approx(-0.5));
  REQUIRE(testTMI(p1) == Approx(-0.5));
  REQUIRE(testTMI(p2) == Approx(-0.5));
  REQUIRE(testTMI(p3) == Approx(0.5));
  // 1x3x1
  elemVals = std::vector<double>(3, 1.);
  nodalVals = std::vector<double>(16, 1.);
  testTMI = VoxelInterp(nodalVals, 1, 3, 1, 1., 1., 1., 0.5);
  REQUIRE(testTMI(p0) == Approx(-0.5));
  REQUIRE(testTMI(p1) == Approx(-0.5));
  REQUIRE(testTMI(p2) == Approx(-0.5));
  REQUIRE(testTMI(p3) == Approx(0.5));
  // 2x2x1
  elemVals = std::vector<double>(4, 1.);
  nodalVals = std::vector<double>(18, 1.);
  testTMI = VoxelInterp(nodalVals, 2, 2, 1, 1., 1., 1., 0.5);
  REQUIRE(testTMI(p0) == Approx(-0.5));
  REQUIRE(testTMI(p1) == Approx(-0.5));
  REQUIRE(testTMI(p2) == Approx(-0.5));
  REQUIRE(testTMI(p3) == Approx(0.50));
  // 2x2x2
  elemVals = std::vector<double>(8, 1.);
  nodalVals = std::vector<double>(27, 1.);
  testTMI = VoxelInterp(nodalVals, 2, 2, 2, 1., 1., 1., 0.5);
  REQUIRE(testTMI(p0) == Approx(-0.5));
  REQUIRE(testTMI(p1) == Approx(-0.5));
  REQUIRE(testTMI(p2) == Approx(-0.5));
  REQUIRE(testTMI(p3) == Approx(0.5));
  // 1x1x2
  elemVals = std::vector<double>(2, 1.);
  nodalVals = {0., 0., 0., 0., 0.5, 0.5, 0.5, 0.5, 1., 1., 1., 1.};
  testTMI = VoxelInterp(nodalVals, 1, 1, 2, 1., 1., 2., 0.5);
  Point_3_base q0(0.5, 0.5, 0.1), q1(0.5, 0.5, 1.), q2(0.5, 0.5, 1.9);
  REQUIRE(testTMI(q0) == Approx( 0.45));
  REQUIRE(testTMI(q1) == Approx(0.));
  REQUIRE(testTMI(q2) == Approx(-0.45));
}
