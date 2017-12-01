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
#include <vector>

using namespace Topologies;

TEST_CASE("Testing wrapper functions for PugiXML","[PugiXMLWrapper]")
{
	using namespace InputLoader;
	pugi::xml_document xmldoc;
	std::string fileName("testpugixml.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	std::vector<std::string> path(1,"rootnode");
	SECTION("PC data, full path")
	{
		// Test full path
		path.push_back("doubledata");
		REQUIRE(readDoublePCData(xmldoc, path) == Approx(10.));
		path.pop_back();
		path.push_back("unsigneddata");
		REQUIRE(readUnsignedPCData(xmldoc, path) == 3);
		path.pop_back();
		path.push_back("intdata");
		REQUIRE(readIntPCData(xmldoc, path) == -2);
		path.pop_back();
		path.push_back("stringdata");
		REQUIRE(readStringPCData(xmldoc, path) == "This is a string.");
		path.pop_back();
		path.push_back("booldata1");
		REQUIRE(readBoolPCData(xmldoc, path) == true);
		path.pop_back();
		path.push_back("booldata2");
		REQUIRE(readBoolPCData(xmldoc, path) == true);
		path.pop_back();
		path.push_back("booldata3");
		REQUIRE(readBoolPCData(xmldoc, path) == false);
		path.pop_back();
		path.push_back("booldata4");
		REQUIRE(readBoolPCData(xmldoc, path) == false);
		path.pop_back();
		path.push_back("file_name1");
		REQUIRE(readAndCheckFileNamePCData(xmldoc, path) == "testtoo.xml");
		path.pop_back();
		path.push_back("file_name2");
		REQUIRE_THROWS_AS(readAndCheckFileNamePCData(xmldoc, path), ParseException); // File doesn't exist
		// Make sure it's the right exception type
		try{readAndCheckFileNamePCData(xmldoc, path);}
		catch(ParseException pe){REQUIRE(pe.pet == petNonexistentFile);}
		// Check nonexistant node
		path.pop_back();
		path.push_back("nothing");
		REQUIRE_THROWS_AS(readDoublePCData(xmldoc, path), ParseException);
		REQUIRE_THROWS_AS(readUnsignedPCData(xmldoc, path), ParseException);
		REQUIRE_THROWS_AS(readIntPCData(xmldoc, path), ParseException);
		REQUIRE_THROWS_AS(readStringPCData(xmldoc, path), ParseException);
		REQUIRE_THROWS_AS(readBoolPCData(xmldoc, path), ParseException);
		REQUIRE_THROWS_AS(readAndCheckFileNamePCData(xmldoc, path), ParseException);
	}
	// Test 1 node path
	pugi::xml_node rootNode = xmldoc.child("rootnode");
	REQUIRE(rootNode);
	SECTION("PC data, 1 node path")
	{
		REQUIRE(readDoublePCData(rootNode, "doubledata") == Approx(10.));
		REQUIRE(readUnsignedPCData(rootNode, "unsigneddata") == 3);
		REQUIRE(readIntPCData(rootNode, "intdata") == -2);
		REQUIRE(readStringPCData(rootNode, "stringdata") == "This is a string.");
		REQUIRE(readBoolPCData(rootNode, "booldata1") == true);
		REQUIRE(readBoolPCData(rootNode, "booldata2") == true);
		REQUIRE(readBoolPCData(rootNode, "booldata3") == false);
		REQUIRE(readBoolPCData(rootNode, "booldata4") == false);
		REQUIRE(readAndCheckFileNamePCData(rootNode, "file_name1") == "testtoo.xml");
		REQUIRE_THROWS_AS(readAndCheckFileNamePCData(rootNode, "file_name2"), ParseException); // File doesn't exist
		// Make sure it's the right exception type
		try{readAndCheckFileNamePCData(rootNode, "file_name2");}
		catch(ParseException pe){REQUIRE(pe.pet == petNonexistentFile);}
		// Check nonexistant node
		REQUIRE_THROWS_AS(readDoublePCData(rootNode, "nothing"), ParseException);
		REQUIRE_THROWS_AS(readUnsignedPCData(rootNode, "nothing"), ParseException);
		REQUIRE_THROWS_AS(readIntPCData(rootNode, "nothing"), ParseException);
		REQUIRE_THROWS_AS(readStringPCData(rootNode, "nothing"), ParseException);
		REQUIRE_THROWS_AS(readBoolPCData(rootNode, "nothing"), ParseException);
		REQUIRE_THROWS_AS(readAndCheckFileNamePCData(rootNode, "nothing"), ParseException);
	}
	SECTION("PC data, vectors")
	{
		std::vector<double> pcDoubleVec = readDoubleVecPCData(rootNode, "doublevec");
		REQUIRE(pcDoubleVec.size() == 3);
		REQUIRE(pcDoubleVec[0] == Approx(1.));
		REQUIRE(pcDoubleVec[1] == Approx(2.));
		REQUIRE(pcDoubleVec[2] == Approx(3.));
		std::vector<unsigned> pcUnsignedVec = readUnsignedVecPCData(rootNode, "unsignedvec");
		REQUIRE(pcUnsignedVec.size() == 4);
		REQUIRE(pcUnsignedVec[0] == 4);
		REQUIRE(pcUnsignedVec[1] == 5);
		REQUIRE(pcUnsignedVec[2] == 6);
		REQUIRE(pcUnsignedVec[3] == 7);
		std::vector<int> pcIntVec = readIntVecPCData(rootNode, "intvec");
		REQUIRE(pcIntVec.size() == 3);
		REQUIRE(pcIntVec[0] == -1);
		REQUIRE(pcIntVec[1] == -2);
		REQUIRE(pcIntVec[2] ==  3);
		// Exceptions
		REQUIRE_THROWS_AS(readDoubleVecPCData(rootNode, "nothing"), ParseException);
		REQUIRE_THROWS_AS(readUnsignedVecPCData(rootNode, "nothing"), ParseException);
		REQUIRE_THROWS_AS(readIntVecPCData(rootNode, "nothing"), ParseException);
	}
	// Test direct node access
	SECTION("PC data, direct node")
	{
		pugi::xml_node pcdnode = rootNode.child("doubledata");
		REQUIRE(readDoublePCData(pcdnode) == Approx(10.));
		pcdnode = rootNode.child("unsigneddata");
		REQUIRE(readUnsignedPCData(pcdnode) == 3);
		pcdnode = rootNode.child("intdata");
		REQUIRE(readIntPCData(pcdnode) == -2);
		pcdnode = rootNode.child("stringdata");
		REQUIRE(readStringPCData(pcdnode) == "This is a string.");
		pcdnode = rootNode.child("booldata1");
		REQUIRE(readBoolPCData(pcdnode) == true);
		pcdnode = rootNode.child("booldata2");
		REQUIRE(readBoolPCData(pcdnode) == true);
		pcdnode = rootNode.child("booldata3");
		REQUIRE(readBoolPCData(pcdnode) == false);
		pcdnode = rootNode.child("booldata4");
		REQUIRE(readBoolPCData(pcdnode) == false);
		pcdnode = rootNode.child("file_name1");
		REQUIRE(readAndCheckFileNamePCData(pcdnode) == "testtoo.xml");
		pcdnode = rootNode.child("file_name2");
		REQUIRE_THROWS_AS(readAndCheckFileNamePCData(pcdnode), ParseException); // File doesn't exist
		// Make sure it's the right exception type
		try{readAndCheckFileNamePCData(pcdnode);}
		catch(ParseException pe){REQUIRE(pe.pet == petNonexistentFile);}
		// Check nonexistant node
		pcdnode = rootNode.child("nothing");
		REQUIRE_FALSE(pcdnode);
		REQUIRE_THROWS_AS(readDoublePCData(pcdnode), ParseException);
		REQUIRE_THROWS_AS(readUnsignedPCData(pcdnode), ParseException);
		REQUIRE_THROWS_AS(readIntPCData(pcdnode), ParseException);
		REQUIRE_THROWS_AS(readStringPCData(pcdnode), ParseException);
		REQUIRE_THROWS_AS(readBoolPCData(pcdnode), ParseException);
	}
	SECTION("Attributes")
	{
		std::string attrname = "attr";
		path.push_back("doubleattr");
		REQUIRE(readDoubleAttribute(xmldoc, path, attrname) == Approx(10.));
		path.pop_back();
		path.push_back("unsignedattr");
		REQUIRE(readUnsignedAttribute(xmldoc, path, attrname) == 3);
		path.pop_back();
		path.push_back("intattr");
		REQUIRE(readIntAttribute(xmldoc, path, attrname) == -2);
		path.pop_back();
		path.push_back("stringattr");
		REQUIRE(readStringAttribute(xmldoc, path, attrname) == "This is a string.");
		path.pop_back();
		path.push_back("boolattr");
		REQUIRE(readBoolAttribute(xmldoc, path, "attr1") == true);
		REQUIRE(readBoolAttribute(xmldoc, path, "attr2") == true);
		REQUIRE(readBoolAttribute(xmldoc, path, "attr3") == false);
		REQUIRE(readBoolAttribute(xmldoc, path, "attr4") == false);
		path.pop_back();
		path.push_back("filenameattr");
		REQUIRE(readAndCheckFileNameAttribute(xmldoc, path, "exists") == "testga.xml");
		REQUIRE_THROWS_AS(readAndCheckFileNameAttribute(xmldoc, path, "doesntexist"), ParseException);
		// Make sure it's the right exception type
		try{readAndCheckFileNameAttribute(xmldoc, path, "doesntexist");}
		catch(ParseException pe){REQUIRE(pe.pet == petNonexistentFile);}
		// Check nonexistant attribute
		attrname = "nothing";
		REQUIRE_THROWS_AS(readDoubleAttribute(xmldoc, path, attrname), ParseException);
		REQUIRE_THROWS_AS(readUnsignedAttribute(xmldoc, path, attrname), ParseException);
		REQUIRE_THROWS_AS(readIntAttribute(xmldoc, path, attrname), ParseException);
		REQUIRE_THROWS_AS(readStringAttribute(xmldoc, path, attrname), ParseException);
		REQUIRE_THROWS_AS(readBoolAttribute(xmldoc, path, attrname), ParseException);
		REQUIRE_THROWS_AS(readAndCheckFileNameAttribute(xmldoc, path, attrname), ParseException);
	}
	SECTION("Attributes, 1 node path")
	{
		std::string attrname = "attr";
		REQUIRE(readDoubleAttribute(rootNode, "doubleattr", attrname) == Approx(10.));
		REQUIRE(readUnsignedAttribute(rootNode, "unsignedattr", attrname) == 3);
		REQUIRE(readIntAttribute(rootNode, "intattr", attrname) == -2);
		REQUIRE(readStringAttribute(rootNode, "stringattr", attrname) == "This is a string.");
		REQUIRE(readBoolAttribute(rootNode, "boolattr", "attr1") == true);
		REQUIRE(readBoolAttribute(rootNode, "boolattr", "attr2") == true);
		REQUIRE(readBoolAttribute(rootNode, "boolattr", "attr3") == false);
		REQUIRE(readBoolAttribute(rootNode, "boolattr", "attr4") == false);
		REQUIRE(readAndCheckFileNameAttribute(rootNode, "filenameattr", "exists") == "testga.xml");
		REQUIRE_THROWS_AS(readAndCheckFileNameAttribute(rootNode, "filenameattr", "doesntexist"), ParseException);
		// Make sure it's the right exception type
		try{readAndCheckFileNameAttribute(rootNode, "filenameattr", "doesntexist");}
		catch(ParseException pe){REQUIRE(pe.pet == petNonexistentFile);}
		// Check nonexistant node
		REQUIRE_THROWS_AS(readDoubleAttribute(rootNode, "doubleattr", "nothing"), ParseException);
		REQUIRE_THROWS_AS(readUnsignedAttribute(rootNode, "doubleattr", "nothing"), ParseException);
		REQUIRE_THROWS_AS(readIntAttribute(rootNode, "doubleattr", "nothing"), ParseException);
		REQUIRE_THROWS_AS(readStringAttribute(rootNode, "doubleattr", "nothing"), ParseException);
		REQUIRE_THROWS_AS(readBoolAttribute(rootNode, "doubleattr", "nothing"), ParseException);
		REQUIRE_THROWS_AS(readAndCheckFileNameAttribute(xmldoc, "doubleattr", "nothing"), ParseException);
	}
	SECTION("Attributes, direct node")
	{
		std::string attrname = "attr";
		pugi::xml_node pcdnode = rootNode.child("doubleattr");
		REQUIRE(readDoubleAttribute(pcdnode, attrname) == Approx(10.));
		pcdnode = rootNode.child("unsignedattr");
		REQUIRE(readUnsignedAttribute(pcdnode, attrname) == 3);
		pcdnode = rootNode.child("intattr");
		REQUIRE(readIntAttribute(pcdnode, attrname) == -2);
		pcdnode = rootNode.child("stringattr");
		REQUIRE(readStringAttribute(pcdnode, attrname) == "This is a string.");
		pcdnode = rootNode.child("boolattr");
		REQUIRE(readBoolAttribute(pcdnode, "attr1") == true);
		REQUIRE(readBoolAttribute(pcdnode, "attr2") == true);
		REQUIRE(readBoolAttribute(pcdnode, "attr3") == false);
		REQUIRE(readBoolAttribute(pcdnode, "attr4") == false);
		pcdnode = rootNode.child("filenameattr");
		REQUIRE(readAndCheckFileNameAttribute(pcdnode, "exists") == "testga.xml");
		REQUIRE_THROWS_AS(readAndCheckFileNameAttribute(pcdnode, "doesntexist"), ParseException);
		// Make sure it's the right exception type
		try{readAndCheckFileNameAttribute(pcdnode, "doesntexist");}
		catch(ParseException pe){REQUIRE(pe.pet == petNonexistentFile);}
		// Check nonexistant node
		pcdnode = rootNode.child("doubleattr");
		attrname = "nothing";
		REQUIRE(pcdnode);
		REQUIRE_THROWS_AS(readDoubleAttribute(pcdnode, attrname), ParseException);
		REQUIRE_THROWS_AS(readUnsignedAttribute(pcdnode, attrname), ParseException);
		REQUIRE_THROWS_AS(readIntAttribute(pcdnode, attrname), ParseException);
		REQUIRE_THROWS_AS(readStringAttribute(pcdnode, attrname), ParseException);
		REQUIRE_THROWS_AS(readBoolAttribute(pcdnode, attrname), ParseException);
		REQUIRE_THROWS_AS(readAndCheckFileNameAttribute(pcdnode, attrname), ParseException);
	}
	SECTION("Node path")
	{
		pugi::xml_node pathNode = rootNode.child("GET").child("Node").child("path");
		REQUIRE(pathNode);
		std::vector<std::string> pathVec = getNodePath(pathNode);
		REQUIRE(pathVec.size() == 4);
		REQUIRE(pathVec[0] == "rootnode");
		REQUIRE(pathVec[1] == "GET");
		REQUIRE(pathVec[2] == "Node");
		REQUIRE(pathVec[3] == "path");
	}
}

