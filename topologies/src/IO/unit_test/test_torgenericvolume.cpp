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

#include "inputloaderrep.h"
#include "catch.hpp"
#include <string>

using namespace Topologies;

TEST_CASE("Testing input parsing using TORGenericVolume","[TORGenericVolume]")
{
	using namespace InputLoader;
	RepNodeInfo testRNI;
	REQUIRE(testRNI.getNodeName() == "representation");
	pugi::xml_document xmldoc;
	std::string fileName("testpix.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("representation");
	REQUIRE(rootNode);
	TORGenericVolume defaultParser(testRNI.getTypeName());
	SECTION("Pixel, no default values")
	{
		testRNI.parse(rootNode, "testpix.xml");
		TORGenericVolume testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getDimensions(0) == Approx(10.));
		REQUIRE(testParser.getDimensions(1) == Approx(1.));
		REQUIRE(testParser.getDiscSizes(0) == 300);
		REQUIRE(testParser.getDiscSizes(1) == 30);
		REQUIRE(testParser.getMET() == metQuad);
		REQUIRE(testParser.getMinDensity() == Approx(10.));
		REQUIRE(testParser.getThreshold() == Approx(50.));
		REQUIRE(testParser.getPenalPower() == Approx(30.));
		REQUIRE(testParser.getFiltRad() == defaultParser.getFiltRad());
		REQUIRE(testParser.getBetaHeavi() == defaultParser.getBetaHeavi());
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Heaviside2D, with default values")
	{
		testRNI.parse(rootNode, "testpix.xml");
		TORGenericVolume testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getDimensions(0) == Approx(1.));
		REQUIRE(testParser.getDimensions(1) == Approx(1.));
		REQUIRE(testParser.getDiscSizes(0) == 30);
		REQUIRE(testParser.getDiscSizes(1) == 30);
		REQUIRE(testParser.getMET() == metTri);
		REQUIRE(testParser.getMinDensity() == defaultParser.getMinDensity());
		REQUIRE(testParser.getThreshold() == defaultParser.getThreshold());
		REQUIRE(testParser.getPenalPower() == defaultParser.getPenalPower());
		REQUIRE(testParser.getFiltRad() == Approx(1.1));
		REQUIRE(testParser.getBetaHeavi() == Approx(0.1));
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Voxel, no defaults")
	{
		testRNI.parse(rootNode, "testpix.xml");
		TORGenericVolume testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getDimensions(0) == Approx(1.));
		REQUIRE(testParser.getDimensions(1) == Approx(0.2));
		REQUIRE(testParser.getDimensions(2) == Approx(1.));
		REQUIRE(testParser.getDiscSizes(0) == 20);
		REQUIRE(testParser.getDiscSizes(1) == 4);
		REQUIRE(testParser.getDiscSizes(2) == 20);
		REQUIRE(testParser.getMET() == metHex);
		REQUIRE(testParser.getMinDensity() == Approx(1.));
		REQUIRE(testParser.getThreshold() == Approx(1.));
		REQUIRE(testParser.getPenalPower() == Approx(2.));
		REQUIRE(testParser.getFiltRad() == defaultParser.getFiltRad());
		REQUIRE(testParser.getBetaHeavi() == defaultParser.getBetaHeavi());
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Heaviside3D, with default values")
	{
		testRNI.parse(rootNode, "testpix.xml");
		TORGenericVolume testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getDimensions(0) == Approx(1.));
		REQUIRE(testParser.getDimensions(1) == Approx(2.));
		REQUIRE(testParser.getDimensions(2) == Approx(3.));
		REQUIRE(testParser.getDiscSizes(0) == 2);
		REQUIRE(testParser.getDiscSizes(1) == 4);
		REQUIRE(testParser.getDiscSizes(2) == 6);
		REQUIRE(testParser.getMET() == metTet);
		REQUIRE(testParser.getMinDensity() == defaultParser.getMinDensity());
		REQUIRE(testParser.getThreshold() == defaultParser.getThreshold());
		REQUIRE(testParser.getPenalPower() == defaultParser.getPenalPower());
		REQUIRE(testParser.getFiltRad() == defaultParser.getFiltRad());
		REQUIRE(testParser.getBetaHeavi() == defaultParser.getBetaHeavi());
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
	}
}
