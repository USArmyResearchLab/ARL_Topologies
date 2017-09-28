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

TEST_CASE("Testing input parsing using TOOGA","[TOOGA]")
{
	using namespace InputLoader;
	OptNodeInfo testONI;
	REQUIRE(testONI.getNodeName() == "optimizer");
	pugi::xml_document xmldoc;
	std::string fileName("testga.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("optimizer");
	REQUIRE(rootNode);
	TOOGA defaultParser(testONI.getTypeName());
	SECTION("GA, no default values")
	{
		testONI.parse(rootNode, "testga.xml");
		TOOGA testParser(testONI.getTypeName());
		testParser.parseNode(testONI);
		REQUIRE(testParser.getConstraintPenalty() == Approx(9.));
		REQUIRE(testParser.getPenaltyPower() == Approx(10.));
		REQUIRE(testParser.getMutationRange() == Approx(5.));
		REQUIRE(testParser.getCrossoverRate() == Approx(3.));
		REQUIRE(testParser.getMutationRate() == Approx(4.));
		REQUIRE(testParser.getPopSize() == 2);
		REQUIRE(testParser.getNumGens() == 1);
		REQUIRE(testParser.getNumTourn() == 6);
		REQUIRE(testParser.getNumElite() == 7);
		REQUIRE(testParser.getMutationRadius() == 8);
		REQUIRE(testParser.getNumGoals() == 1);
		REQUIRE(testParser.getSharingRadius() == Approx(11.));
		REQUIRE_FALSE(testParser.getIsPareto());
		REQUIRE(testParser.getGoalWeights().empty());
	}
	rootNode = rootNode.next_sibling("optimizer");
	REQUIRE(rootNode);
	SECTION("PGA, with default values")
	{
		testONI.parse(rootNode, "testga.xml");
		TOOGA testParser(testONI.getTypeName());
		testParser.parseNode(testONI);
		REQUIRE(testParser.getConstraintPenalty() == Approx(defaultParser.getConstraintPenalty()));
		REQUIRE(testParser.getPenaltyPower() == Approx(defaultParser.getPenaltyPower()));
		REQUIRE(testParser.getMutationRange() == Approx(16.));
		REQUIRE(testParser.getCrossoverRate() == Approx(14.));
		REQUIRE(testParser.getMutationRate() == Approx(15.));
		REQUIRE(testParser.getPopSize() == 13);
		REQUIRE(testParser.getNumGens() == 12);
		REQUIRE(testParser.getNumTourn() == defaultParser.getNumTourn());
		REQUIRE(testParser.getNumElite() == defaultParser.getNumElite());
		REQUIRE(testParser.getMutationRadius() == defaultParser.getMutationRadius());
		REQUIRE(testParser.getNumGoals() == 3);
		REQUIRE(testParser.getSharingRadius() == defaultParser.getSharingRadius());
		REQUIRE(testParser.getIsPareto());
		const std::vector<double>& goalWeights = testParser.getGoalWeights();
		REQUIRE(goalWeights.size() == 3);
		REQUIRE(goalWeights[0] == 1);
		REQUIRE(goalWeights[1] == 2);
		REQUIRE(goalWeights[2] == 3);
	}
}
