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

#include "exotxtmeshloader.h"
#include "tomesh.h"
#include "catch.hpp"
#include <string>

using namespace Topologies;

typedef std::vector<std::pair<unsigned, std::vector<std::size_t>>> NodeSet;
typedef std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> SideSet;

TEST_CASE("Testing 2d, tri element functions in ExoTxtMeshLoader namespace","[ExoTxtMeshLoader]")
{
	using namespace InputLoader::ExoTxtMeshLoader;
	const std::string fileName("testtri.txt");
	REQUIRE(readMeshDimension(fileName) == 2);
	std::unique_ptr<TOMesh2D> testMesh = loadExoTxtFileTri(fileName);
	REQUIRE(testMesh);
	SECTION("Test mesh")
	{
		REQUIRE(testMesh->getNumNodes() == 41);
		REQUIRE(testMesh->getNumElements() == 60);
	}
	SECTION("Test node set")
	{
		NodeSet testns = loadNodeSets(fileName);
		REQUIRE(testns.size() == 1); // 1 node set
		REQUIRE(testns[0].first == 1); // id of 1
		REQUIRE(testns[0].second.size() == 6); // 6 nodes in node set 1
		// Check location, should all be on x=-0.5 plane
		for(auto it = testns[0].second.begin(); it != testns[0].second.end(); ++it)
			REQUIRE(testMesh->getNode2D(*it).x() == Approx(-0.5));
	}
	SECTION("Test side set")
	{
		SideSet testss = loadSideSets(fileName, testMesh.get());
		REQUIRE(testss.size() == 1); // 1 side set
		REQUIRE(testss[0].first == 1); // id of 1
		REQUIRE(testss[0].second.size() == 5); // 5 edges in side set 1
		for(auto it = testss[0].second.begin(); it != testss[0].second.end(); ++it)
			REQUIRE(it->size() == 2);
		// Check location, should all be on x=0.5 line
		for(auto it1 = testss[0].second.begin(); it1 != testss[0].second.end(); ++it1)
			for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
				REQUIRE(testMesh->getNode2D(*it2).x() == Approx(0.5));
	}
}

TEST_CASE("Testing 2d, quad element functions in ExoTxtMeshLoader namespace","[ExoTxtMeshLoader]")
{
	using namespace InputLoader::ExoTxtMeshLoader;
	const std::string fileName("testquad.txt");
	REQUIRE(readMeshDimension(fileName) == 2);
	std::unique_ptr<TOMesh2D> testMesh = loadExoTxtFileTri(fileName);
	REQUIRE(testMesh);
	SECTION("Test mesh")
	{
		REQUIRE(testMesh->getNumNodes() == 36);
		REQUIRE(testMesh->getNumElements() == 25);
	}
	SECTION("Test node set")
	{
		NodeSet testns = loadNodeSets(fileName);
		REQUIRE(testns.size() == 1); // 1 node set
		REQUIRE(testns[0].first == 1); // id of 1
		REQUIRE(testns[0].second.size() == 6); // 6 nodes in node set 1
		// Check location, should all be on x=-0.5 plane
		for(auto it = testns[0].second.begin(); it != testns[0].second.end(); ++it)
			REQUIRE(testMesh->getNode2D(*it).x() == Approx(-0.5));
	}
	SECTION("Test side set")
	{
		SideSet testss = loadSideSets(fileName, testMesh.get());
		REQUIRE(testss.size() == 1); // 1 side set
		REQUIRE(testss[0].first == 1); // id of 1
		REQUIRE(testss[0].second.size() == 5); // 5 edges in side set 1
		for(auto it = testss[0].second.begin(); it != testss[0].second.end(); ++it)
			REQUIRE(it->size() == 2);
		// Check location, should all be on x=0.5 line
		for(auto it1 = testss[0].second.begin(); it1 != testss[0].second.end(); ++it1)
			for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
				REQUIRE(testMesh->getNode2D(*it2).x() == Approx(0.5));
	}
}

TEST_CASE("Testing two blocks in 2D","[ExoTxtMeshLoader]")
{
	using namespace InputLoader::ExoTxtMeshLoader;
	const std::string fileName("test2blocks.txt");
	REQUIRE(readMeshDimension(fileName) == 2);
	std::unique_ptr<TOMesh2D> testMesh = loadExoTxtFileTri(fileName);
	REQUIRE(testMesh);
	SECTION("Test mesh")
	{
		REQUIRE(testMesh->getNumNodes() == 240);
		REQUIRE(testMesh->getNumElements() == 210);
	}
	SECTION("Test block ids")
	{
		for(std::size_t k = 0; k < 105; ++k)
			REQUIRE(testMesh->getMatID(k) == 1);
		for(std::size_t k = 105; k < 210; ++k)
			REQUIRE(testMesh->getMatID(k) == 2);
	}
}

TEST_CASE("Testing 3d, tet functions in ExoTxtMeshLoader namespace","[ExoTxtMeshLoader]")
{
	using namespace InputLoader::ExoTxtMeshLoader;
	const std::string fileName("testtet.txt");
	REQUIRE(readMeshDimension(fileName) == 3);
	std::unique_ptr<TOMesh3D> testMesh = loadExoTxtFileTet(fileName);
	REQUIRE(testMesh);
	SECTION("Test mesh")
	{
		REQUIRE(testMesh->getNumNodes() == 296);
		REQUIRE(testMesh->getNumElements() == 1175);
	}
	SECTION("Test node set")
	{
		NodeSet testns = loadNodeSets(fileName);
		REQUIRE(testns.size() == 1); // 1 node set
		REQUIRE(testns[0].first == 1); // id of 1
		REQUIRE(testns[0].second.size() == 39); // 39 nodes in node set 1
		// Check location, should all be on x=-0.5 plane
		for(auto it = testns[0].second.begin(); it != testns[0].second.end(); ++it)
			REQUIRE(testMesh->getNode3D(*it).x() == Approx(-0.5));
	}
	SECTION("Test side set")
	{
		SideSet testss = loadSideSets(fileName, testMesh.get());
		REQUIRE(testss.size() == 1); // 1 side set
		REQUIRE(testss[0].first == 1); // id of 1
		REQUIRE(testss[0].second.size() == 56); // 56 faces in side set 1
		for(auto it = testss[0].second.begin(); it != testss[0].second.end(); ++it)
			REQUIRE(it->size() == 3);
		// Check location, should all be on x=0.5 plane
		for(auto it1 = testss[0].second.begin(); it1 != testss[0].second.end(); ++it1)
			for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
				REQUIRE(testMesh->getNode3D(*it2).x() == Approx(0.5));
	}
}

TEST_CASE("Testing 3d, hex element functions in ExoTxtMeshLoader namespace","[ExoTxtMeshLoader]")
{
	using namespace InputLoader::ExoTxtMeshLoader;
	const std::string fileName("testhex.txt");
	REQUIRE(readMeshDimension(fileName) == 3);
	std::unique_ptr<TOMesh3D> testMesh = loadExoTxtFileTet(fileName);
	REQUIRE(testMesh);
	SECTION("Test mesh")
	{
		REQUIRE(testMesh->getNumNodes() == 216);
		REQUIRE(testMesh->getNumElements() == 125);
	}
	SECTION("Test node set")
	{
		NodeSet testns = loadNodeSets(fileName);
		REQUIRE(testns.size() == 1); // 1 node set
		REQUIRE(testns[0].first == 1); // id of 1
		REQUIRE(testns[0].second.size() == 36); // 36 nodes in node set 1
		// Check location, should all be on x=-0.5 plane
		for(auto it = testns[0].second.begin(); it != testns[0].second.end(); ++it)
			REQUIRE(testMesh->getNode3D(*it).x() == Approx(-0.5));
	}
	SECTION("Test side set")
	{
		SideSet testss = loadSideSets(fileName, testMesh.get());
		REQUIRE(testss.size() == 1); // 1 side set
		REQUIRE(testss[0].first == 1); // id of 1
		REQUIRE(testss[0].second.size() == 25); // 25 faces in side set 1
		for(auto it = testss[0].second.begin(); it != testss[0].second.end(); ++it)
			REQUIRE(it->size() == 4);
		// Check location, should all be on x=0.5 plane
		for(auto it1 = testss[0].second.begin(); it1 != testss[0].second.end(); ++it1)
			for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
				REQUIRE(testMesh->getNode3D(*it2).x() == Approx(0.5));
	}
}

