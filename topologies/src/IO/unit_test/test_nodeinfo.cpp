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

#include "inputloader.h"
#include "catch.hpp"
#include <string>

using namespace Topologies;

TEST_CASE("Testing RepNodeInfo class","[RepNodeInfo]")
{
	using namespace InputLoader;
	RepNodeInfo testRNI;
	REQUIRE(testRNI.getNodeName() == "representation");
	pugi::xml_document xmldoc;
	std::string fileName("nodeinfo.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("representation");
	REQUIRE(rootNode);
	SECTION("Pixel")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "pixel");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testpix.xml");
		REQUIRE(testRNI.getType() == tortPixel);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Heaviside2D")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "heaviside2d");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "nodeinfo.xml");
		REQUIRE(testRNI.getType() == tortHeaviside2D);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Voxel")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "voxel");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testpix.xml");
		REQUIRE(testRNI.getType() == tortVoxel);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Heaviside3D")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "heaviside3d");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testpix.xml");
		REQUIRE(testRNI.getType() == tortHeaviside3D);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("CSG2D")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "csg2d");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testcsg.xml");
		REQUIRE(testRNI.getType() == tortCSG2D);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("CSG3D")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "csg3d");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testcsg.xml");
		REQUIRE(testRNI.getType() == tortCSG3D);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Mesh2D")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "mesh2d");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testmesh2d.xml");
		REQUIRE(testRNI.getType() == tortMesh2D);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Mesh3D")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "mesh3d");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testmesh2d.xml");
		REQUIRE(testRNI.getType() == tortMesh3D);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("HeavisideMesh2D")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "heavisidemesh2d");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testmesh2d.xml");
		REQUIRE(testRNI.getType() == tortHeavisideMesh2D);
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("HeavisideMesh3D")
	{
		testRNI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testRNI.getTypeName() == "heavisidemesh3d");
		REQUIRE(testRNI.getTag() == "testrep");
		std::vector<std::string> path = testRNI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "representation");
		REQUIRE(testRNI.getFileName() == "testmesh2d.xml");
		REQUIRE(testRNI.getType() == tortHeavisideMesh3D);
	}
}

TEST_CASE("Testing OptNodeInfo classes","[RepNodeInfo]")
{
	using namespace InputLoader;
	OptNodeInfo testONI;
	REQUIRE(testONI.getNodeName() == "optimizer");
	pugi::xml_document xmldoc;
	std::string fileName("nodeinfo.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("optimizer");
	REQUIRE(rootNode);
	SECTION("OC")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "oc");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testtoo.xml");
		REQUIRE(testONI.getType() == tootOC);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("MMA")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "mma");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "nodeinfo.xml");
		REQUIRE(testONI.getType() == tootMMA);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("BFGS")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "bfgs");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testtoo.xml");
		REQUIRE(testONI.getType() == tootBFGS);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("GA")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "ga");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testga.xml");
		REQUIRE(testONI.getType() == tootGA);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("PGA")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "pga");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testga.xml");
		REQUIRE(testONI.getType() == tootPGA);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("GD")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "gd");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testtoo.xml");
		REQUIRE(testONI.getType() == tootGD);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("Chain")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "chain");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testchain.xml");
		REQUIRE(testONI.getType() == tootChain);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("Refine")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "refine");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testtoo.xml");
		REQUIRE(testONI.getType() == tootRefine);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("Convert")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "convert");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testtoo.xml");
		REQUIRE(testONI.getType() == tootConvert);
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("Continuation")
	{
		testONI.parse(rootNode, "nodeinfo.xml");
		REQUIRE(testONI.getTypeName() == "continuation");
		REQUIRE(testONI.getTag() == "testopt");
		std::vector<std::string> path = testONI.getPath();
		REQUIRE(path.size() == 1);
		REQUIRE(path[0] == "optimizer");
		REQUIRE(testONI.getFileName() == "testchain.xml");
		REQUIRE(testONI.getType() == tootContinuation);
	}
}
