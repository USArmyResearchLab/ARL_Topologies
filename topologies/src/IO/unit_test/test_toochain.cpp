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

TEST_CASE("Testing input parsing using TOOChain","[TOOChain]")
{
	using namespace InputLoader;
	OptNodeInfo testONI;
	REQUIRE(testONI.getNodeName() == "optimizer");
	pugi::xml_document xmldoc;
	std::string fileName("testchain.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("optimizer");
	REQUIRE(rootNode);
	testONI.parse(rootNode, fileName);
	TOOChain testParser;
	testParser.parseNode(testONI);
	// Check each entry
	auto it = testParser.oniBegin(); 
	REQUIRE(it != testParser.oniEnd());
	REQUIRE(it->getTypeName() == "oc");
	REQUIRE(it->getType() == tootOC);
	++it;
	REQUIRE(it != testParser.oniEnd());
	REQUIRE(it->getTypeName() == "refine");
	REQUIRE(it->getType() == tootRefine);
	++it;
	REQUIRE(it != testParser.oniEnd());
	REQUIRE(it->getTypeName() == "convert");
	REQUIRE(it->getType() == tootConvert);
	++it;
	REQUIRE(it != testParser.oniEnd());
	REQUIRE(it->getTypeName() == "continuation");
	REQUIRE(it->getType() == tootContinuation);
	++it;
	REQUIRE(it == testParser.oniEnd());
}
