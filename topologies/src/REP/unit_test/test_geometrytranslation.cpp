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

#include "geometrytranslation.h"
#include "meshtestns.h"
#include "catch.hpp"
#include <algorithm>
#include <iostream>

using namespace Topologies;

void printMesh(const CDT_2& mesh)
{
	std::cout << "Vertices:" << std::endl;
	for(auto it = mesh.points_begin(); it != mesh.points_end(); ++it)
		std::cout << *it << std::endl;
	std::cout << "Faces: " << std::endl;
	for(auto it = mesh.finite_faces_begin(); it != mesh.finite_faces_end(); ++it)
		std::cout << it->vertex(0)->info() << " " << it->vertex(1)->info() << " " << it->vertex(2)->info() << " " << std::endl;
}

TEST_CASE("Testing GeometryTranslation namespace, 2d mesh functions", "[GeometryTranslation]")
{
	using namespace GeometryTranslation;
	using namespace MeshTestNS;
	MesherData md; // Use default values
	std::vector<GenericMaterial> matVec1 = {GenericMaterial({0., 1.})};
	std::vector<GenericMaterial> matVec2 = {GenericMaterial({0., 1.}), GenericMaterial({0., 0.})};
	std::vector<Mesh_Segment_2> segVec = {Mesh_Segment_2(Point_2_base(0., 0.), Point_2_base(1., 0.)),
		Mesh_Segment_2(Point_2_base(1., 0.), Point_2_base(1., 2.)),
		Mesh_Segment_2(Point_2_base(1., 2.), Point_2_base(0., 2.)),
		Mesh_Segment_2(Point_2_base(0., 2.), Point_2_base(0., 0.))};
	// First mesh, no material properties
	CDT_2 testMesh = mesh2D(segVec, md);
	REQUIRE(testMesh.is_valid());
	std::pair<double,double> xrange = getXRange(testMesh);
	std::pair<double,double> yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(1.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(2.));
	REQUIRE(getMeshArea(testMesh) == Approx(2.));
	// With matvec
	testMesh = mesh2D(segVec, matVec1, md);
	REQUIRE(testMesh.is_valid());
	xrange = getXRange(testMesh);
	yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(1.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(2.));
	REQUIRE(getMeshArea(testMesh) == Approx(2.));
	REQUIRE(computeMaterialPropertyTotal(testMesh, 0, 2) == 0.);
	REQUIRE(computeMaterialPropertyTotal(testMesh, 1, 2) == 2.);
	// With matvec holding 2 properties
	testMesh = mesh2D(segVec, matVec2, md);
	REQUIRE(testMesh.is_valid());
	xrange = getXRange(testMesh);
	yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(1.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(2.));
	REQUIRE(getMeshArea(testMesh) == Approx(2.));
	REQUIRE(computeMaterialPropertyTotal(testMesh, 0, 2) == 0.);
	REQUIRE(computeMaterialPropertyTotal(testMesh, 1, 2) == 2.);
	// Two polygons
	std::vector<Mesh_Segment_2> segVec2 = {Mesh_Segment_2(Point_2_base(0.4, 0.4), Point_2_base(0.6, 0.4)),
		Mesh_Segment_2(Point_2_base(0.6, 0.4), Point_2_base(0.6, 0.6)),
		Mesh_Segment_2(Point_2_base(0.6, 0.6), Point_2_base(0.4, 0.6)),
		Mesh_Segment_2(Point_2_base(0.4, 0.6), Point_2_base(0.4, 0.4))};
	std::vector<std::vector<Mesh_Segment_2>> segVV = {segVec, segVec2};
	// No material properties
	testMesh = mesh2D(segVV, md);
	REQUIRE(testMesh.is_valid());
	xrange = getXRange(testMesh);
	yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(1.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(2.));
	REQUIRE(getMeshArea(testMesh) == Approx(1.96));
	// One material property
	testMesh = mesh2D(segVec, matVec1, md);
	REQUIRE(testMesh.is_valid());
	xrange = getXRange(testMesh);
	yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(1.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(2.));
	REQUIRE(getMeshArea(testMesh) == Approx(2.));
	REQUIRE(computeMaterialPropertyTotal(testMesh, 0, 2) == Approx(0.));
	REQUIRE(computeMaterialPropertyTotal(testMesh, 1, 2) == Approx(2.));
	// 2 material properties (should be the same as above)
	testMesh = mesh2D(segVV, matVec2, md);
	REQUIRE(testMesh.is_valid());
	xrange = getXRange(testMesh);
	yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(1.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(2.));
	REQUIRE(getMeshArea(testMesh) == Approx(1.96));
	REQUIRE(computeMaterialPropertyTotal(testMesh, 0, 2) == Approx(0.));
	REQUIRE(computeMaterialPropertyTotal(testMesh, 1, 2) == Approx(1.96));
	// pixel mesh functions
	// First, no material properties
	testMesh = mesh2Dpixel(3, 2, 1., 2., md);
	REQUIRE(testMesh.is_valid());
	xrange = getXRange(testMesh);
	yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(1.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(2.));
	REQUIRE(getMeshArea(testMesh) == Approx(2.));
	// With material vector
	std::vector<GenericMaterial> matVec3 = {GenericMaterial({1.}), GenericMaterial({2.}), 
		GenericMaterial({3.}), GenericMaterial({4.}), GenericMaterial({5.}), GenericMaterial({6.})};
	testMesh = mesh2Dpixel(2, 3, 2., 1., matVec3, md);
	REQUIRE(testMesh.is_valid());
	xrange = getXRange(testMesh);
	yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(2.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(1.));
	REQUIRE(getMeshArea(testMesh) == Approx(2.));
	REQUIRE(computeMaterialPropertyTotal(testMesh, 0, 1) == Approx(7.));
	// with opt val vector
	std::vector<double> optValVec = {0., 1., 2., 3., 4., 5.};
	testMesh = mesh2Dpixel(3, 2, 1., 2., optValVec, md);
	REQUIRE(testMesh.is_valid());
	xrange = getXRange(testMesh);
	yrange = getYRange(testMesh);
	REQUIRE(xrange.first == Approx(0.));
	REQUIRE(xrange.second == Approx(1.));
	REQUIRE(yrange.first == Approx(0.));
	REQUIRE(yrange.second == Approx(2.));
	REQUIRE(getMeshArea(testMesh) == Approx(2.));
	REQUIRE(computeOptValTotal(testMesh) == Approx(5.));
}

TEST_CASE("Testing GeometryTranslation namespace, 3d mesh functions", "[GeometryTranslation]")
{
	using namespace GeometryTranslation;
	using namespace MeshTestNS;
	// Build 1x1x1 cube out of triangles
	std::vector<Point_3_base> ptVec({Point_3_base(0., 0., 0.), Point_3_base(1., 0., 0.),
																		Point_3_base(0., 1., 0.), Point_3_base(1., 1., 0.),
																		Point_3_base(0., 0., 1.), Point_3_base(1., 0., 1.),
																		Point_3_base(0., 1., 1.), Point_3_base(1., 1., 1.)});
	std::vector<std::vector<std::size_t>> connectivity({{0, 2, 1}, {2, 3, 1}, {4, 5, 6}, {6, 5, 7},
																											{0, 5, 4}, {0, 1, 5}, {2, 6, 7}, {2, 7, 3},
																											{1, 7, 5}, {1, 3, 7}, {0, 4, 6}, {0, 6, 2}});
	TOMesh3DSurface tom3s(ptVec, connectivity);
	// Create a Mesh_polyhedron_3
	Mesh_polyhedron_3 MP3;
	GeometryTranslation::PolyhedronBuilderFromTOMesh builder(tom3s);
	MP3.delegate(builder);
	// Mesh
	MesherData md; // Use default values
	C3t3 mesh = mesh3D(MP3, md);
	// Get point ranges
	std::vector<double> bbox = getBoundBox(mesh);
	REQUIRE(bbox[0] == Approx(0.));
	REQUIRE(bbox[1] == Approx(1.));
	REQUIRE(bbox[2] == Approx(0.));
	REQUIRE(bbox[3] == Approx(1.));
	REQUIRE(bbox[4] == Approx(0.));
	REQUIRE(bbox[5] == Approx(1.));
	// Check volume
	REQUIRE(getMeshVolume(mesh) == Approx(1.));
}

TEST_CASE("Boundary mesh functions", "[GeometryTranslation]")
{
	using namespace GeometryTranslation;
	using namespace MeshTestNS;
	std::vector<Mesh_Segment_2> segVec = {Mesh_Segment_2(Point_2_base(-0.5, -1.), Point_2_base(0.5, -1.)),
		Mesh_Segment_2(Point_2_base(0.5, -1.), Point_2_base(0.5, 1.)),
		Mesh_Segment_2(Point_2_base(0.5, 1.), Point_2_base(-0.5, 1.)),
		Mesh_Segment_2(Point_2_base(-0.5, 1.), Point_2_base(-0.5, -1.))};
	std::vector<Mesh_Segment_2> segVec2 = {Mesh_Segment_2(Point_2_base(0., 0.), Point_2_base(0.1, 0.)),
		Mesh_Segment_2(Point_2_base(0.1, 0.), Point_2_base(0.1, 0.1)),
		Mesh_Segment_2(Point_2_base(0.1, 0.1), Point_2_base(0., 0.1)),
		Mesh_Segment_2(Point_2_base(0., 0.1), Point_2_base(0., 0.))};
	std::vector<std::vector<Mesh_Segment_2>> segVV;
	REQUIRE(computeSignedArea(segVec) == Approx(2.));
	SECTION("Segments from points")
	{
		// Test mesh segment formation
		std::vector<Point_2_base> ptVec = {Point_2_base(-0.5, -1.), Point_2_base(0.5, -1.),
			Point_2_base(0.5, 1.), Point_2_base(-0.5, 1.)};
		segVec.clear();
		addPointsToSegments(ptVec, segVec);
		REQUIRE(segVec.size() == 4);
		REQUIRE(computeSignedArea(segVec) == Approx(2.));
		addPointsToSegments(ptVec, segVV);
		REQUIRE(segVV.size() == 1);
		REQUIRE(segVV[0].size() == 4);
		REQUIRE(computeSignedArea(segVV[0]) == Approx(2.));
		// Similar function, but wipes segVec first
		generateSegment2FromPoints(ptVec, segVec);
		REQUIRE(segVec.size() == 4);
		REQUIRE(computeSignedArea(segVec) == Approx(2.));
		// Get points from segments
		generatePointsFromSegment2(segVec, ptVec);
		REQUIRE(ptVec.size() == 4);
		generateSegment2FromPoints(ptVec, segVec);
		REQUIRE(segVec.size() == 4);
		REQUIRE(computeSignedArea(segVec) == Approx(2.));
	}
	SECTION("Segment ordering")
	{
		// Test all permutations of segVec
		std::vector<std::size_t> kvec = {0, 1, 2, 3};
		do
		{
			std::vector<Mesh_Segment_2> curVec = {segVec[kvec[0]], segVec[kvec[1]], segVec[kvec[2]], segVec[kvec[3]]};
			orderMeshSegments(curVec, segVV);
			REQUIRE(segVV.size() == 1);
			REQUIRE(segVV[0].size() == 4);
			REQUIRE(computeSignedArea(segVV[0]) == Approx(2.));
		} while(std::next_permutation(kvec.begin(), kvec.end()));
		// Add a hole polygon
		segVec.insert(segVec.end(), segVec2.begin(), segVec2.end());
		kvec = {0, 1, 2, 3, 4, 5, 6, 7};
		do
    {
      std::vector<Mesh_Segment_2> curVec(8);
			for(std::size_t k = 0; k < segVec.size(); ++k)
				curVec[k] = segVec[kvec[k]];
      orderMeshSegments(curVec, segVV);
      REQUIRE(segVV.size() == 2);
      REQUIRE(segVV[0].size() == 4);
			REQUIRE(segVV[1].size() == 4);
			// The order could be different, so need to check first
			double area0 = computeSignedArea(segVV[0]), area1 = computeSignedArea(segVV[1]);
			if(area0 > area1)
				std::swap(area0, area1);
			REQUIRE(area0 == Approx(-0.01));
			REQUIRE(area1 == Approx(2.));
			// Test flatten
			flattenMeshSegments(segVV, curVec);
			REQUIRE(curVec.size() == 8);
			
		} while(std::next_permutation(kvec.begin(), kvec.end()));
		// Test flatten
		std::vector<Mesh_Segment_2> curVec;
		flattenMeshSegments(segVV, curVec);
		REQUIRE(curVec.size() == 8);
		// Check boundary validity
		REQUIRE(checkValidityAsBoundaryMesh(segVV[0]));
		REQUIRE(checkValidityAsBoundaryMesh(segVV[1]));
		REQUIRE(checkValidityAsBoundaryMesh(segVec2));
		segVec2.push_back(Mesh_Segment_2(Point_2_base(0., 0.), Point_2_base(-0.1, 0.)));
		REQUIRE(!checkValidityAsBoundaryMesh(segVec2));
	}
	SECTION("CGAL polygon type")
	{
		// Segs to polygon
		CGAL::Polygon_2<Mesh_K> poly;
		meshSegs2Polygon(segVec, poly);
		REQUIRE(poly.area() == Approx(2.));
		// Polygon to segs
		polygon2MeshSegs(poly, segVec);
		REQUIRE(computeSignedArea(segVec) == Approx(2.));
		// Nested polygon
		segVec.insert(segVec.end(), segVec2.begin(), segVec2.end());
		std::vector<CGAL::Polygon_2<Mesh_K>> polyVec;
		meshSegs2PolygonWithHoles(segVec, poly, polyVec);
		REQUIRE(polyVec.size() == 1);
		REQUIRE(poly.area() == Approx(2.));
		REQUIRE(polyVec[0].area() == Approx(-0.01));
		// Point projection
		segVec.erase(segVec.begin() + 4, segVec.end());
		Point_2_base testP(0., 2.);
		projectPointToPolygon(testP, poly, 1e-14);
		REQUIRE(testP.x() == Approx(0.));
		REQUIRE(testP.y() == Approx(1.));
		testP = Point_2_base(1., 0.);
		projectPointToPolygon(testP, poly, 1e-14);
		REQUIRE(testP.x() == Approx(0.5));
		REQUIRE(testP.y() == Approx(0.0));
		testP = Point_2_base(-1., -2.);
		projectPointToPolygon(testP, poly, 1e-14);
		REQUIRE(testP.x() == Approx(-0.5));
		REQUIRE(testP.y() == Approx(-1.0));
		// Component vector
		std::vector<double> compVec = {-0.5, -1., 0.5, -1., 0.5, 1., -0.5, 1.};
		componentVec2Polygon(compVec, poly);
		REQUIRE(poly.area() == Approx(2.));
		std::vector<double> compVec2;
		polygon2ComponentVec(poly, compVec2);
		REQUIRE(compVec == compVec2);
		// Partition
		segVec.erase(segVec.end()-1);
		segVec.push_back(Mesh_Segment_2(Point_2_base(-0.5, 1.), Point_2_base(0., 0.)));
		segVec.push_back(Mesh_Segment_2(Point_2_base(0., 0.), Point_2_base(-0.5, -1.)));
		std::vector<std::vector<Mesh_Segment_2> > segVV;
		partitionMeshSegPolys(segVec, segVV);
		double area = 0.;
		REQUIRE(segVV.size() == 2);
		REQUIRE((computeSignedArea(segVV[0]) + computeSignedArea(segVV[1])) == Approx(1.5));
	}
}

TEST_CASE("Testing GeometryTranslation namespace, extrusion", "[GeometryTranslation]")
{
	using namespace GeometryTranslation;
	using namespace MeshTestNS;
	MesherData md; // Use default values
	std::vector<Mesh_Segment_2> segVec = {Mesh_Segment_2(Point_2_base(0., 0.), Point_2_base(1., 0.)),
		Mesh_Segment_2(Point_2_base(1., 0.), Point_2_base(1., 2.)),
		Mesh_Segment_2(Point_2_base(1., 2.), Point_2_base(0., 2.)),
		Mesh_Segment_2(Point_2_base(0., 2.), Point_2_base(0., 0.))};
	// First mesh, no material properties
	CDT_2 testMesh = mesh2D(segVec, md);
	REQUIRE(testMesh.is_valid());
	TOMesh3DSurface tom3s = linearExtrude(testMesh, 0.1);
	// Check number of points & triangles
	REQUIRE(2*testMesh.number_of_vertices() == tom3s.getNumNodes());
	REQUIRE(2*testMesh.number_of_faces() < tom3s.getNumElements());
	// Check surface area
	REQUIRE(computeSurfaceArea(tom3s) == Approx(4.6));
}

