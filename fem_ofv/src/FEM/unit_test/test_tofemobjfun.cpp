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

#include "tofemobjfun.h"
#include "topoptrep.h"
#include "torfactory.h"
#include "inputloaderrep.h"
#include "catch.hpp"
#include <memory>

using namespace Topologies;

InputLoader::RepNodeInfo loadRNI(const std::string& fileName)
{
	// Loads an XML document and parses the first <representation> node seen
	pugi::xml_document xmldoc;
	xmldoc.load_file(fileName.c_str());
	assert(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("representation");
	InputLoader::RepNodeInfo testRNI;
	testRNI.parse(rootNode, fileName);
	return testRNI;
}

void testGradient(TopOptRep& testTOR, TOFEMObjFun& testObjFun)
{
	std::pair<std::vector<double>, bool> res1, res2;
	testObjFun.g(testTOR, res1);
	REQUIRE(res1.second);
	testObjFun.TopOptObjFun::g(testTOR, res2); //Calls finite difference version
	REQUIRE(res2.second);
	REQUIRE(res1.first.size() == res2.first.size());
	for(std::size_t k = 0; k < res1.first.size(); ++k)
		REQUIRE(res1.first[k] == Approx(res2.first[k]).epsilon(1e-4));
}

void testConstraintGradient(TopOptRep& testTOR, TOFEMObjFun& testObjFun)
{
	std::pair<std::vector<double>, bool> res1, res2;
	testObjFun.gc(testTOR, res1);
	REQUIRE(res1.second);
	testObjFun.TopOptObjFun::gc(testTOR, res2); //Calls finite difference version
	REQUIRE(res2.second);
	REQUIRE(res1.first.size() == res2.first.size());
	for(std::size_t k = 0; k < res1.first.size(); ++k)
		REQUIRE(res1.first[k] == Approx(res2.first[k]).epsilon(1e-4));
}

TEST_CASE("Testing TOFEMObjFun with VolMesh2D", "[TOFEMObjFun]")
{
	// Load objective function
	TOFEMObjFun testObjFun("tof_exo_2.xml");
	// Load TORep
	InputLoader::RepNodeInfo testRNI = loadRNI("testmesh2d.xml");
	std::unique_ptr<TopOptRep> upVM2D = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upVM2D);
	upVM2D->initialize(0.5);
	// Test gradient
	testGradient(*upVM2D, testObjFun);
	testConstraintGradient(*upVM2D, testObjFun);
}

TEST_CASE("Testing TOFEMObjFun with VolMesh2D and GMSH input", "[TOFEMObjFun]")
{
	// Load objective function
	TOFEMObjFun testObjFun("tof_gmsh_2.xml");
	// Load TORep
	InputLoader::RepNodeInfo testRNI = loadRNI("testmesh2dGMSH.xml");
	std::unique_ptr<TopOptRep> upVM2D = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upVM2D);
	upVM2D->initialize(0.5);
	// Test gradient
	testGradient(*upVM2D, testObjFun);
	testConstraintGradient(*upVM2D, testObjFun);
}

TEST_CASE("Testing TOFEMObjFun HeavisideMesh2D", "[TOFEMObjFun]")
{
	// Load objective function
  TOFEMObjFun testObjFun("tof_exo_heavi2.xml");
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testheavimesh2d.xml");
	std::unique_ptr<TopOptRep> upHeavi = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upHeavi);
	upHeavi->initialize(0.5);
	// Test gradient
	testGradient(*upHeavi, testObjFun);
	testConstraintGradient(*upHeavi, testObjFun);
}

TEST_CASE("Testing TOFEMObjFun with VolMesh3D", "[TOFEMObjFun]")
{
	// Load objective function
  TOFEMObjFun testObjFun("tof_exo_3.xml");
	InputLoader::RepNodeInfo testRNI = loadRNI("testmesh3d.xml");
	std::unique_ptr<TopOptRep> upVM3D = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upVM3D);
	upVM3D->initialize(0.5);
	// Test gradient
	testGradient(*upVM3D, testObjFun);
	testConstraintGradient(*upVM3D, testObjFun);
}

TEST_CASE("Testing TOFEMObjFun HeavisideMesh3D", "[TOFEMObjFun]")
{
	// Load objective function
	TOFEMObjFun testObjFun("tof_exo_heavi3.xml");
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testheavimesh3d.xml");
	std::unique_ptr<TopOptRep> upHeavi = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upHeavi);
	upHeavi->initialize(0.5);
	// Test gradient
	testGradient(*upHeavi, testObjFun);
	testConstraintGradient(*upHeavi, testObjFun);
}

