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

TEST_CASE("Testing input parsing using TORCSGTree","[TORCSGTree]")
{
	using namespace InputLoader;
	RepNodeInfo testRNI;
	REQUIRE(testRNI.getNodeName() == "representation");
	pugi::xml_document xmldoc;
	std::string fileName("testcsg.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("representation");
	REQUIRE(rootNode);
	TORCSGTree defaultParser(testRNI.getTypeName());
	SECTION("Pixel, no default values")
	{
		testRNI.parse(rootNode, "testcsg.xml");
		TORCSGTree testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getRegionDimensions(0) == Approx(1.));
		REQUIRE(testParser.getRegionDimensions(1) == Approx(2.));
		REQUIRE(testParser.getMeshSizes(0) == 3);
		REQUIRE(testParser.getMeshSizes(1) == 4);
		REQUIRE(testParser.getShapeNums(0) == 5);
		REQUIRE(testParser.getShapeNums(1) == 6);
		REQUIRE(testParser.getNumPointsPerShape() == 10);
		REQUIRE(testParser.getMinDensity() == Approx(1.));
		REQUIRE(testParser.getUseAffine() == true);
		REQUIRE(testParser.getShapesAreHoles() == true);
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Heaviside2D, with default values")
	{
		testRNI.parse(rootNode, "testcsg.xml");
		TORCSGTree testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getRegionDimensions(0) == Approx(1.));
    REQUIRE(testParser.getRegionDimensions(1) == Approx(2.));
		REQUIRE(testParser.getRegionDimensions(2) == Approx(3.));
    REQUIRE(testParser.getMeshSizes(0) == 4);
    REQUIRE(testParser.getMeshSizes(1) == 5);
		REQUIRE(testParser.getMeshSizes(2) == 6);
    REQUIRE(testParser.getShapeNums(0) == 7);
    REQUIRE(testParser.getShapeNums(1) == 8);
		REQUIRE(testParser.getShapeNums(2) == 9);
    REQUIRE(testParser.getNumPointsPerShape() == 10);
    REQUIRE(testParser.getMinDensity() == defaultParser.getMinDensity());
    REQUIRE(testParser.getUseAffine() == defaultParser.getUseAffine());
    REQUIRE(testParser.getShapesAreHoles() == defaultParser.getShapesAreHoles());
    REQUIRE(testParser.getRepName() == testRNI.getTypeName());
	}
}
