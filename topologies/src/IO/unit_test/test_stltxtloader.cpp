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

#include "stltxtloader.h"
#include "tomesh.h"
#include "catch.hpp"
#include <string>
#include <CGAL/Triangle_3.h>

using namespace Topologies;

TEST_CASE("Testing STL load","[STLTxtLoader]")
{
	std::string fileName("cube.stl");
	SECTION("Test CGAL Polyhedron_3")
	{
		Mesh_polyhedron_3 mp3 = InputLoader::STLTxtLoader::loadFileCGAL(fileName);
		// Check surface area
		double sa = 0.;
		for(auto fit = mp3.facets_begin(); fit != mp3.facets_end(); ++fit)
		{
			Mesh_K::Triangle_3 t3(fit->halfedge()->vertex()->point(), fit->halfedge()->next()->vertex()->point(),
          fit->halfedge()->opposite()->vertex()->point());
			sa += sqrt(t3.squared_area());
		}
		REQUIRE(sa == Approx(600.));
	}
	SECTION("Test TOMesh")
	{
		std::unique_ptr<TOMesh3DSurface> upTOM = InputLoader::STLTxtLoader::loadFileTOMesh(fileName);
		// Check surface area
		double sa = 0.;
		for(std::size_t k = 0; k < upTOM->getNumElements(); ++k)
		{
			const std::vector<std::size_t>& connVec = upTOM->getElementConnectivity(k);
			REQUIRE(connVec.size() == 3);
			Mesh_K::Triangle_3 t3(upTOM->getNode3D(connVec[0]), upTOM->getNode3D(connVec[1]), upTOM->getNode3D(connVec[2]));
			sa += sqrt(t3.squared_area());
		}
		REQUIRE(sa == Approx(600.));
	}
}

