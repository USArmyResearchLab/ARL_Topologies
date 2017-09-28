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

#include "inputloadertopopt.h"
#include "catch.hpp"
#include <string>

using namespace Topologies;

TEST_CASE("Testing main input parsing","[Topologies]")
{
	using namespace InputLoader;
	pugi::xml_document xmldoc;
	std::string fileName("testtopopt.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	TopologiesParser testParser;
	SECTION("Test 1")
	{
		testParser.parse(xmldoc);
		REQUIRE(testParser.getTOFFileName() == "tof_geo.xml");
		REQUIRE(testParser.getTOFSharedLibName() == "libfemofv.so");
		REQUIRE(testParser.getNumProcessorsPerOFV() == 10);
		const RepNodeInfo& testRNI = testParser.getRepNodeInfo();
		REQUIRE(testRNI.getType() == tortPixel);
		REQUIRE(testRNI.getTag() == "rep1");
		const OptNodeInfo& testONI = testParser.getOptNodeInfo();
		REQUIRE(testONI.getType() == tootOC);
		REQUIRE(testONI.getTag() == "opt1");
		const InitialGuess& testIG = testParser.getInitialGuess();
		REQUIRE(testParser.getInitialGuess().getInitialGuessType() == igtFile);
		std::size_t numOutputs = std::distance(testParser.outputBegin(), testParser.outputEnd());
		REQUIRE(numOutputs == 2);
	}
	pugi::xml_node rootNode = xmldoc.child("topologies_test2");
	TopologiesParser defaultParser;
	SECTION("Test 2")
	{
		testParser.parse(rootNode);
		REQUIRE(testParser.getNumProcessorsPerOFV() == defaultParser.getNumProcessorsPerOFV());
		std::size_t numOutputs = std::distance(testParser.outputBegin(), testParser.outputEnd());
		REQUIRE(numOutputs == 0);
		const InitialGuess& testIG = testParser.getInitialGuess();
		InitialGuess defaultIG;
		REQUIRE(testIG.getInitialGuessType() == defaultIG.getInitialGuessType());
		REQUIRE(testIG.getConstant() == defaultIG.getConstant());
	}
}

TEST_CASE("Testing initial guess parsing","[InitialGuess]")
{
	using namespace InputLoader;
	pugi::xml_document xmldoc;
	std::string fileName("testinitialguess.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	InitialGuess testParser, defaultParser;
	SECTION("Test constant")
	{
		testParser.parse(xmldoc);
		REQUIRE(testParser.getInitialGuessType() == igtConstant);
		REQUIRE(testParser.getConstant() == Approx(1.));
	}
	pugi::xml_node rootNode = xmldoc.child("initial_guess2");
	SECTION("Test constant with noise")
	{
		testParser.parse(rootNode);
		REQUIRE(testParser.getInitialGuessType() == igtConstantWithNoise);
		REQUIRE(testParser.getConstant() == Approx(2.));
		REQUIRE(testParser.getRandRange().first == Approx(3.));
		REQUIRE(testParser.getRandRange().second == Approx(4.));
	}
	rootNode = xmldoc.child("initial_guess3");
	SECTION("Test random")
	{
		testParser.parse(rootNode);
		REQUIRE(testParser.getInitialGuessType() == igtRandom);
		REQUIRE(testParser.getRandRange().first == Approx(5.));
		REQUIRE(testParser.getRandRange().second == Approx(6.));
	}
	rootNode = xmldoc.child("initial_guess4");
	SECTION("Test file")
	{
		testParser.parse(rootNode);
		REQUIRE(testParser.getInitialGuessType() == igtFile);
		REQUIRE(testParser.getFileName() == "testfile.txt");
	}
}

TEST_CASE("Testing output parsing","[Output]")
{
	using namespace InputLoader;
	pugi::xml_document xmldoc;
	std::string fileName("testoutput.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	Output defaultParser;
	SECTION("Volume")
	{
		Output testParser;
		testParser.parse(xmldoc);
		REQUIRE(testParser.getOutputType() == otVolume);
		REQUIRE(testParser.getFileName() == "vol");
		REQUIRE(testParser.getFileFormat() == offVTK);
		REQUIRE(testParser.getExtrusionLength() == Approx(defaultParser.getExtrusionLength()));
		REQUIRE(testParser.getOutputPeriod() == 5);
		REQUIRE(testParser.getOverwrite() == false);
		REQUIRE(testParser.getOutputFinal() == defaultParser.getOutputFinal());
		REQUIRE(testParser.getOutputPeriodic() == true);
	}
	pugi::xml_node rootNode = xmldoc.child("output2");
	SECTION("Surface")
	{
		Output testParser;
		testParser.parse(rootNode);
		REQUIRE(testParser.getOutputType() == otSurface);
		REQUIRE(testParser.getFileName() == "segs");
		REQUIRE(testParser.getFileFormat() == offMatlab);
		REQUIRE(testParser.getOutputPeriod() == defaultParser.getOutputPeriod());
		REQUIRE(testParser.getOverwrite() == defaultParser.getOverwrite());
		REQUIRE(testParser.getOutputFinal() == true);
		REQUIRE(testParser.getOutputPeriodic() == defaultParser.getOutputPeriodic());
	}
	rootNode = xmldoc.child("output3");
	SECTION("Raw data")
	{
		Output testParser;
		testParser.parse(rootNode);
		REQUIRE(testParser.getOutputType() == otRawData);
		REQUIRE(testParser.getFileName() == "data");
		REQUIRE(testParser.getFileFormat() == offDefault);
		REQUIRE(testParser.getOutputPeriod() == 1);
		REQUIRE(testParser.getOverwrite() == true);
		REQUIRE(testParser.getOutputFinal() == false);
		REQUIRE(testParser.getOutputPeriodic() == true);
	}
	rootNode = xmldoc.child("output4");
	SECTION("Extrude")
	{
		Output testParser;
		testParser.parse(rootNode);
		REQUIRE(testParser.getOutputType() == otExtrude);
		REQUIRE(testParser.getFileName() == "extrusion");
		REQUIRE(testParser.getFileFormat() == offSTL);
		REQUIRE(testParser.getExtrusionLength() == Approx(0.1));
		REQUIRE(testParser.getOutputPeriod() == defaultParser.getOutputPeriod());
		REQUIRE(testParser.getOverwrite() == defaultParser.getOverwrite());
		REQUIRE(testParser.getOutputFinal() == defaultParser.getOutputFinal());
		REQUIRE(testParser.getOutputPeriodic() == defaultParser.getOutputPeriodic());
	}
	rootNode = xmldoc.child("output5");
	SECTION("Objective function result")
	{
		Output testParser;
		testParser.parse(rootNode);
		REQUIRE(testParser.getOutputType() == otObjFunRes);
		REQUIRE(testParser.getFileName() == "ofvres");
		REQUIRE(testParser.getFileFormat() == offDefault);
		REQUIRE(testParser.getExtrusionLength() == Approx(defaultParser.getExtrusionLength()));
		REQUIRE(testParser.getOutputPeriod() == defaultParser.getOutputPeriod());
		REQUIRE(testParser.getOverwrite() == defaultParser.getOverwrite());
		REQUIRE(testParser.getOutputFinal() == defaultParser.getOutputFinal());
		REQUIRE(testParser.getOutputPeriodic() == defaultParser.getOutputPeriodic());
	}
}
