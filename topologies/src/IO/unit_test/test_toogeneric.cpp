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

#include "inputloaderopt.h"
#include "catch.hpp"
#include <string>

using namespace Topologies;

TEST_CASE("Testing input parsing using TOOGeneric","[TOOGeneric]")
{
	using namespace InputLoader;
	OptNodeInfo testONI;
	REQUIRE(testONI.getNodeName() == "optimizer");
	pugi::xml_document xmldoc;
	std::string fileName("testtoo.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("optimizer");
	REQUIRE(rootNode);
	TOOGeneric defaultParser(testONI.getTypeName());
	SECTION("OC, no default values")
	{
		testONI.parse(rootNode, "testtoo.xml");
		TOOGeneric testParser(testONI.getTypeName());
		testParser.parseNode(testONI);
		REQUIRE(testParser.getFilterSize() == Approx(1.));
		REQUIRE(testParser.getConstraintPenalty() == defaultParser.getConstraintPenalty());
		REQUIRE(testParser.getPenaltyPower() == defaultParser.getPenaltyPower());
		REQUIRE(testParser.getStepSize() == defaultParser.getStepSize());
		REQUIRE(testParser.getMaxStep() == defaultParser.getMaxStep());
		REQUIRE(testParser.getStopTol() == Approx(2.));
		REQUIRE(testParser.getMaxIters() == 3);
		REQUIRE(testParser.getOptimizerName() == testONI.getTypeName());
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("MMA, with default values")
	{
		testONI.parse(rootNode, "testtoo.xml");
		TOOGeneric testParser(testONI.getTypeName());
		testParser.parseNode(testONI);
		REQUIRE(testParser.getFilterSize() == defaultParser.getFilterSize());
		REQUIRE(testParser.getConstraintPenalty() == defaultParser.getConstraintPenalty());
		REQUIRE(testParser.getPenaltyPower() == defaultParser.getPenaltyPower());
		REQUIRE(testParser.getStepSize() == defaultParser.getStepSize());
		REQUIRE(testParser.getMaxStep() == defaultParser.getMaxStep());
		REQUIRE(testParser.getStopTol() == defaultParser.getStopTol());
		REQUIRE(testParser.getMaxIters() == defaultParser.getMaxIters());
		REQUIRE(testParser.getOptimizerName() == testONI.getTypeName());
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("BFGS, no defaults")
	{
		testONI.parse(rootNode, "testtoo.xml");
		TOOGeneric testParser(testONI.getTypeName());
		testParser.parseNode(testONI);
		REQUIRE(testParser.getFilterSize() == Approx(1.));
		REQUIRE(testParser.getConstraintPenalty() == Approx(4.));
		REQUIRE(testParser.getPenaltyPower() == Approx(5.));
		REQUIRE(testParser.getStepSize() == defaultParser.getStepSize());
		REQUIRE(testParser.getMaxStep() == defaultParser.getMaxStep());
		REQUIRE(testParser.getStopTol() == Approx(2.));
		REQUIRE(testParser.getMaxIters() == 3);
		REQUIRE(testParser.getOptimizerName() == testONI.getTypeName());
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("GD, with default values")
	{
		testONI.parse(rootNode, "testtoo.xml");
		TOOGeneric testParser(testONI.getTypeName());
		testParser.parseNode(testONI);
		REQUIRE(testParser.getFilterSize() == defaultParser.getFilterSize());
		REQUIRE(testParser.getConstraintPenalty() == defaultParser.getConstraintPenalty());
		REQUIRE(testParser.getPenaltyPower() == defaultParser.getPenaltyPower());
		REQUIRE(testParser.getStepSize() == Approx(2.));
		REQUIRE(testParser.getMaxStep() == Approx(1.));
		REQUIRE(testParser.getStopTol() == defaultParser.getStopTol());
		REQUIRE(testParser.getMaxIters() == defaultParser.getMaxIters());
		REQUIRE(testParser.getOptimizerName() == testONI.getTypeName());
	}
}
