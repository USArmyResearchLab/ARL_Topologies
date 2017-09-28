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

#include "tomesh.h"
#include "meshtestns.h"
#include "tomeshprocessing.h"
#include "geometrytranslation.h"
#include "catch.hpp"

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

using namespace Topologies;

void addInfoToCDT(CDT_2& inCDT)
{
	// Mark faces
	for(CDT_2::Finite_faces_iterator it = inCDT.finite_faces_begin(); it != inCDT.finite_faces_end(); ++it)
	{
		it->info().nesting_level = -1;
		it->info().optVal = 1.;
	}
	// Add node id to vertices
	unsigned id = 0;
	for(CDT_2::Finite_vertices_iterator it = inCDT.finite_vertices_begin(); it != inCDT.finite_vertices_end(); ++it)
		it->info() = id++;
}

unsigned countFaces(const CDT_2& inCDT)
{
	unsigned nfaces = 0;
	for(CDT_2::Finite_faces_iterator fit = inCDT.finite_faces_begin(); fit != inCDT.finite_faces_end(); ++fit)
	{
		if(fit->is_in_domain())
			++nfaces;
	}
	return nfaces;
}

unsigned countFaces(const C2t3& inMesh)
{
	const SM3_Tr& tr = inMesh.triangulation();
	unsigned nfaces = 0;
	for(SM3_Tr::Finite_facets_iterator fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
	{
		const typename SM3_Tr::Cell_handle cell = fit->first;
		const int& index = fit->second;
		if(cell->is_facet_on_surface(index)==true)
			++nfaces;
	}
	return nfaces;
}

unsigned countVertices(const C2t3& inMesh)
{
	const SM3_Tr& tr = inMesh.triangulation();
	unsigned nverts = 0;
	for(SM3_Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit)
		++nverts;
	return nverts;
}

double sumOptVals(const TOMesh& mesh)
{
	double optSum = 0.;
	for(std::size_t k = 0; k < mesh.getNumElements(); ++k)
		optSum += mesh.getOptVal(k);
	return optSum;
}

// CGAL surface mesh of a sphere (take from examples)
typedef Tr_GT::FT (*Function)(Point_3_base);
typedef CGAL::Implicit_surface_3<Tr_GT, Function> Surface_3;
Tr_GT::FT sphere_function (Point_3_base p) 
{
	const Tr_GT::FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
	return x2 + y2 + z2 - 1.;
}

TEST_CASE("Testing TOMesh in 2D", "[TOMesh]")
{
	using namespace MeshTestNS;
	// Test conversion from CDT_2
	// Generate CDT mesh, uses CGAL example
	CDT_2 testCDT;
	CDT_2::Vertex_handle va = testCDT.insert(Point_2_base(-4, 0));
	CDT_2::Vertex_handle vb = testCDT.insert(Point_2_base( 0,-1));
	CDT_2::Vertex_handle vc = testCDT.insert(Point_2_base( 4, 0));
	CDT_2::Vertex_handle vd = testCDT.insert(Point_2_base( 0, 1));
	testCDT.insert(Point_2_base(2, 0.6)); // Outside region
	testCDT.insert_constraint(va, vb);
	testCDT.insert_constraint(vb, vc);
	testCDT.insert_constraint(vc, vd);
	testCDT.insert_constraint(vd, va);
	Mesher_2 mesher(testCDT);
	mesher.refine_mesh();
	addInfoToCDT(testCDT);
	// Convert to TOMesh
	TOMesh2D testTOM(testCDT);
	REQUIRE(testTOM.dimNum() == 2);
	REQUIRE(testTOM.getNumNodes() == testCDT.number_of_vertices());
	REQUIRE(testTOM.getNumElements() == countFaces(testCDT));
	REQUIRE(getMeshArea(testTOM) == Approx(getMeshArea(testCDT)));
	REQUIRE(sumOptVals(testTOM) == Approx((double)testTOM.getNumElements()));
	// Remesh
	mesher.set_criteria(Criteria(0.125, 0.5));
	mesher.refine_mesh();
	addInfoToCDT(testCDT);
	testTOM = TOMesh2D(testCDT); // Test assignment operator
	REQUIRE(testTOM.dimNum() == 2);
	REQUIRE(testTOM.getNumNodes() == testCDT.number_of_vertices());
	REQUIRE(testTOM.getNumElements() == countFaces(testCDT));
	REQUIRE(getMeshArea(testTOM) == Approx(getMeshArea(testCDT)));
	REQUIRE(sumOptVals(testTOM) == Approx((double)testTOM.getNumElements()));
	// Test other constructors
	// Copy ctor
	TOMesh2D copyTOM(testTOM);
	REQUIRE(copyTOM.dimNum() == 2);
	REQUIRE(copyTOM.getNumNodes() == testCDT.number_of_vertices());
  REQUIRE(copyTOM.getNumElements() == countFaces(testCDT));
	REQUIRE(getMeshArea(copyTOM) == Approx(getMeshArea(testCDT)));
	REQUIRE(sumOptVals(copyTOM) == Approx((double)copyTOM.getNumElements()));
}

TEST_CASE("Testing TOMesh in 3D, volume", "[TOMesh]")
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
	REQUIRE(computeSurfaceArea(tom3s) == Approx(6.));
  // Create a Mesh_polyhedron_3
  Mesh_polyhedron_3 MP3;
  GeometryTranslation::PolyhedronBuilderFromTOMesh builder(tom3s);
  MP3.delegate(builder);
  // Mesh
  MesherData md;
/*	md.tetMeshCellSize = 0.5;
	md.tetMeshEdgeSize = 0.5;
	md.tetMeshFacetAngle = 25.;
	md.tetMeshFacetSize = 0.5;
	md.tetMeshFacetDistance = 0.125;
	md.tetMeshCellRadiusEdgeRatio = 3.;*/
  C3t3 mesh = mesh3D(MP3, md);
	TOMesh3D testTOM(mesh);
	// Test result
	REQUIRE(testTOM.dimNum() == 3);
	REQUIRE(testTOM.getNumElements() == mesh.number_of_cells_in_complex());
	REQUIRE(testTOM.getNumNodes() == mesh.triangulation().number_of_vertices());
	REQUIRE(getMeshArea(testTOM) == Approx(getMeshVolume(mesh)));
	REQUIRE(sumOptVals(testTOM) == Approx(0.));
	testTOM.setOptVals(std::vector<double>(testTOM.getNumElements(), 1.)); // Set optvals
	REQUIRE(sumOptVals(testTOM) == Approx((double)testTOM.getNumElements()));
	// Assignment operator
	testTOM = TOMesh3D(mesh);
	REQUIRE(testTOM.dimNum() == 3);
	REQUIRE(testTOM.getNumElements() == mesh.number_of_cells_in_complex());
	REQUIRE(testTOM.getNumNodes() == mesh.triangulation().number_of_vertices());
	REQUIRE(getMeshArea(testTOM) == Approx(getMeshVolume(mesh)));
	REQUIRE(sumOptVals(testTOM) == Approx(0.));
	testTOM.setOptVals(std::vector<double>(testTOM.getNumElements(), 1.)); // Set optvals
	REQUIRE(sumOptVals(testTOM) == Approx((double)testTOM.getNumElements()));
	// Copy ctor
	TOMesh3D copyTOM(testTOM);
	REQUIRE(copyTOM.dimNum() == 3);
	REQUIRE(copyTOM.getNumElements() == mesh.number_of_cells_in_complex());
	REQUIRE(copyTOM.getNumNodes() == mesh.triangulation().number_of_vertices());
	REQUIRE(getMeshArea(copyTOM) == Approx(getMeshVolume(mesh)));
	REQUIRE(sumOptVals(copyTOM) == Approx((double)testTOM.getNumElements()));
}

TEST_CASE("Testing TOMesh in 3D, surface mesh", "[TOMesh]")
{
	SM3_Tr tr;
	C2t3 c2t3(tr);
	Surface_3 surface(sphere_function, Tr_GT::Sphere_3(CGAL::ORIGIN, 2.));
	CGAL::Surface_mesh_default_criteria_3<SM3_Tr> criteria(30., 0.1, 0.1);
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
	// Test creation
	TOMesh3DSurface testTOM(c2t3);
	REQUIRE(testTOM.dimNum() == 3);
	REQUIRE(testTOM.getNumElements() == countFaces(c2t3));
	REQUIRE(testTOM.getNumNodes() == countVertices(c2t3));
	REQUIRE(sumOptVals(testTOM) == Approx(0.));
	// Test assignment operator
	testTOM = TOMesh3DSurface(c2t3);
	REQUIRE(testTOM.dimNum() == 3);
	REQUIRE(testTOM.getNumElements() == countFaces(c2t3));
	REQUIRE(testTOM.getNumNodes() == countVertices(c2t3));
	REQUIRE(sumOptVals(testTOM) == Approx(0.));
	// Test copy ctor
	TOMesh3DSurface copyTOM(testTOM);
	REQUIRE(copyTOM.dimNum() == 3);
	REQUIRE(copyTOM.getNumElements() == countFaces(c2t3));
	REQUIRE(copyTOM.getNumNodes() == countVertices(c2t3));
	REQUIRE(sumOptVals(copyTOM) == Approx(0.));
}

