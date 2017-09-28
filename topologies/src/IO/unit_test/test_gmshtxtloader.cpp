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

#include "gmshtxtloader.h"
#include "tomesh.h"
#include "catch.hpp"
#include <string>

using namespace Topologies;

typedef std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> PEVector;

TEST_CASE("Testing 2d, tri element functions in GMSHTxtMeshLoader namespace","[GMSHTxtMeshLoader]")
{
	using namespace InputLoader::GMSHTxtMeshLoader;
	const std::string fileName("testtri.msh");
	std::unique_ptr<TOMesh2D> testMesh = loadGMSHTxtFile2D(fileName);
	REQUIRE(testMesh);
	SECTION("Test mesh")
	{
		REQUIRE(testMesh->getNumNodes() == 145);
		REQUIRE(testMesh->getNumElements() == 248);
	}
	SECTION("Test physical entity")
	{
		PEVector testpe = loadBoundaryPhysicalEntity(fileName, 2);
		REQUIRE(testpe.size() == 2); // 2 physical entities
		REQUIRE(testpe[0].first == 7); // id of 7
		REQUIRE(testpe[0].second.size() == 1); // 1 nodes in physical entity 1
		// Check location, should be at (0.5,-0.5)
		REQUIRE(testpe[0].second[0].size() == 1);
		REQUIRE(testMesh->getNode2D(testpe[0].second[0][0]).x() == Approx(0.5));
		REQUIRE(testMesh->getNode2D(testpe[0].second[0][0]).y() == Approx(-0.5));
		// Second group
		REQUIRE(testpe[1].first == 8); // id of 8
		REQUIRE(testpe[1].second.size() == 10); // 10 edges in pe 8
		for(auto it = testpe[1].second.begin(); it != testpe[1].second.end(); ++it)
			REQUIRE(it->size() == 2); // All should be line segments
		// Check location, should all be on x=-0.5 line
		for(auto it1 = testpe[1].second.begin(); it1 != testpe[1].second.end(); ++it1)
			for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
				REQUIRE(testMesh->getNode2D(*it2).x() == Approx(-0.5));
	}
}

TEST_CASE("Testing two blocks in 2D","[GMSHTxtMeshLoader]")
{
	using namespace InputLoader::GMSHTxtMeshLoader;
	const std::string fileName("test2blocks.msh");
	std::unique_ptr<TOMesh2D> testMesh = loadGMSHTxtFile2D(fileName);
	REQUIRE(testMesh);
	SECTION("Test mesh")
	{
		REQUIRE(testMesh->getNumNodes() == 52);
		REQUIRE(testMesh->getNumElements() == 80);
	}
	SECTION("Test block ids")
	{
		for(std::size_t k = 0; k < 40; ++k)
			REQUIRE(testMesh->getMatID(k) == 13);
		for(std::size_t k = 41; k < 80; ++k)
			REQUIRE(testMesh->getMatID(k) == 12);
	}
}

TEST_CASE("Testing 3d, tet functions in GMSHTxtMeshLoader namespace","[GMSHTxtMeshLoader]")
{
	using namespace InputLoader::GMSHTxtMeshLoader;
	const std::string fileName("testtet.msh");
	std::unique_ptr<TOMesh3D> testMesh = loadGMSHTxtFile3D(fileName);
	REQUIRE(testMesh);
	SECTION("Test mesh")
	{
		REQUIRE(testMesh->getNumNodes() == 238);
		REQUIRE(testMesh->getNumElements() == 726);
	}
	SECTION("Test node set")
	{
		PEVector testpe = loadBoundaryPhysicalEntity(fileName, 3);
		REQUIRE(testpe.size() == 3); // 3 pes
		// Check first pe
		REQUIRE(testpe[0].first == 28); // id of 28
		REQUIRE(testpe[0].second.size() == 68); // 68 faces in side set 1
		for(auto it = testpe[0].second.begin(); it != testpe[0].second.end(); ++it)
			REQUIRE(it->size() == 3); // All triangles
		// Check location, should all be on x=-0.5 plane
		for(auto it1 = testpe[0].second.begin(); it1 != testpe[0].second.end(); ++it1)
			for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
				REQUIRE(testMesh->getNode3D(*it2).x() == Approx(-0.5));
		// Check next pe
		REQUIRE(testpe[1].first == 29); // id of 29
		REQUIRE(testpe[1].second.size() == 68); // 68 faces in side set 1
		for(auto it = testpe[1].second.begin(); it != testpe[1].second.end(); ++it)
			REQUIRE(it->size() == 3); // All triangles
		// Check location, should all be on x=0.5 plane
		for(auto it1 = testpe[1].second.begin(); it1 != testpe[1].second.end(); ++it1)
			for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
				REQUIRE(testMesh->getNode3D(*it2).x() == Approx(0.5));
		// Check last PE
		REQUIRE(testpe[2].first == 30); // id of 30
		REQUIRE(testpe[2].second.size() == 1); // 39 nodes in node set 1
		// Check location, should be at (0.5,0.5,-0.5)
		REQUIRE(testpe[2].second[0].size() == 1);
		REQUIRE(testMesh->getNode3D(testpe[2].second[0][0]).x() == Approx(0.5));
		REQUIRE(testMesh->getNode3D(testpe[2].second[0][0]).y() == Approx(0.5));
		REQUIRE(testMesh->getNode3D(testpe[2].second[0][0]).z() == Approx(-0.5));
	}
}

