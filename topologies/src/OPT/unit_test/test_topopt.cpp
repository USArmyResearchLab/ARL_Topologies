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

#include "topoptoc.h"
#include "topoptgd.h"
#include "topoptnlopt.h"
#include "topoptga.h"
#include "topoptchain.h"
#include "topoptfactory.h"
#include "catch.hpp"
#include "objfuntest2.h"
#include "reptest.h"
#include "inputloaderopt.h"
#include <memory>
#include <string>

using namespace Topologies;

namespace Topologies{class OutputHandler;}

InputLoader::OptNodeInfo loadONI(const std::string& fileName)
{
	InputLoader::OptNodeInfo testONI;
	pugi::xml_document xmldoc;
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("optimizer");
	REQUIRE(rootNode);
	testONI.parse(rootNode, fileName);
	return testONI;
}

InputLoader::TOOGeneric loadOptParams(const std::string& fileName)
{
	InputLoader::OptNodeInfo testONI = loadONI(fileName);
	InputLoader::TOOGeneric testParser(testONI.getTypeName());
	testParser.parseNode(testONI);
	return testParser;
}

InputLoader::TOOGA loadOptParamsGA(const std::string& fileName)
{
	InputLoader::OptNodeInfo testONI = loadONI(fileName);
	InputLoader::TOOGA testParser(testONI.getTypeName());
	testParser.parseNode(testONI);
	return testParser;
}

InputLoader::TOOChain loadOptParamsChain(const std::string& fileName)
{
	InputLoader::OptNodeInfo testONI = loadONI(fileName);
	InputLoader::TOOChain testParser;
	testParser.parseNode(testONI);
	return testParser;
}

void testTopOpt(TopOpt& testOpt, const TopOptRep& ig, const std::vector<double>& checkRes, double chkEps = 1e-5)
{
	std::unique_ptr<TopOptRep> upRes = testOpt.optimize(ig);
	std::vector<double> resVec;
	upRes->getRealRep(resVec);
	REQUIRE(resVec.size() == checkRes.size());
	for(std::size_t k = 0; k < resVec.size(); ++k)
		REQUIRE(resVec[k] == Approx(checkRes[k]).epsilon(chkEps));
}

TEST_CASE("Testing gradient descent method","[TopOptGD]")
{
	double target = 0.5;
	TOTestObjFun2 testObjFun(target, 0.);
	TestRep testTORep;
	testTORep.randomize(); // Sets initial guess to random values between 0 & 1
	REQUIRE(testObjFun.gradCheck(testTORep));
	std::vector<OutputHandler*> ohVec;
	TopOptGD testOpt(loadOptParams("gdinput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {target, target});
	// Test 2nd target
	target = 0.75;
	testObjFun = TOTestObjFun2(target, 0.);
	testOpt = TopOptGD(loadOptParams("gdinput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {target, target});
	// Test factory creation
	std::unique_ptr<TopOpt> upTO = TopOptFactory::createTopOpt(loadONI("gdinput.xml"), &testObjFun, ohVec);
	testTopOpt(*upTO, testTORep, {target, target});
}

TEST_CASE("Testing optimality criterion method","[TopOptOC]")
{
	double target = 0.5, constraint = 0.25;
	TOTestObjFun2 testObjFun(target, constraint, false);
	TestRep testTORep;
	testTORep.setRealRep({0.2, 0.2});
	std::vector<OutputHandler*> ohVec;
	TopOptOC testOpt(loadOptParams("ocinput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint/2., constraint/2.}, 1e-4);
	// Test factory creation
  std::unique_ptr<TopOpt> upTO = TopOptFactory::createTopOpt(loadONI("ocinput.xml"), &testObjFun, ohVec);
	testTopOpt(*upTO, testTORep, {constraint/2., constraint/2.}, 1e-4);
}

TEST_CASE("Testing method of moving asymptotes","[TopOptMMA]")
{
	double target = 0.5, constraint = 0.25;
	TOTestObjFun2 testObjFun(target, constraint);
	TestRep testTORep;
	testTORep.setRealRep({0.2, 1.});
	std::vector<OutputHandler*> ohVec;
	TopOptNLOpt testOpt(tootMMA, loadOptParams("mmainput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint, target});
	// Test 2nd target
	target = 0.75;
	testObjFun = TOTestObjFun2(target, constraint);
	testOpt = TopOptNLOpt(tootMMA, loadOptParams("mmainput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint, target});
	// Test factory creation
	std::unique_ptr<TopOpt> upTO = TopOptFactory::createTopOpt(loadONI("mmainput.xml"), &testObjFun, ohVec);
	testTopOpt(*upTO, testTORep, {constraint, target});
}

TEST_CASE("Testing BFGS method","[TopOptBFGS]")
{
	double target = 0.5, constraint = 0.;
	TOTestObjFun2 testObjFun(target, constraint);
	TestRep testTORep;
	testTORep.setRealRep({0.1, 1.});
	std::vector<OutputHandler*> ohVec;
	TopOptNLOpt testOpt(tootBFGS, loadOptParams("bfgsinput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint, target});
	// Test 2nd target
	target = 0.75;
	testObjFun = TOTestObjFun2(target, constraint);
	testOpt = TopOptNLOpt(tootBFGS, loadOptParams("bfgsinput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint, target});
	// Test 3rd
	target = 0.5;
	constraint = 0.25;
	testObjFun = TOTestObjFun2(target, constraint);
	testTORep.setRealRep({0.2, 1.});
	testTopOpt(testOpt, testTORep, {constraint, target});
	// Test factory creation
	std::unique_ptr<TopOpt> upTO = TopOptFactory::createTopOpt(loadONI("bfgsinput.xml"), &testObjFun, ohVec);
	testTopOpt(*upTO, testTORep, {constraint, target});
}

TEST_CASE("Testing GA method","[TopOptGA]")
{
	double target = 0.5, constraint = 0.;
	TOTestObjFun2 testObjFun(target, constraint);
	TestRep testTORep;
	testTORep.setRealRep({0.1, 1.});
	std::vector<OutputHandler*> ohVec;
	TopOptRealGA testOpt(loadOptParamsGA("gainput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint, target}, 1e-3);
	// Test 2nd target
	target = 0.75;
	testObjFun = TOTestObjFun2(target, constraint);
	testOpt = TopOptRealGA(loadOptParamsGA("gainput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint, target}, 1e-3);
	// Test factory creation
	std::unique_ptr<TopOpt> upTO = TopOptFactory::createTopOpt(loadONI("gainput.xml"), &testObjFun, ohVec);
	testTopOpt(*upTO, testTORep, {constraint, target}, 1e-3);
}

TEST_CASE("Testing chain method","[TopOptChain]")
{
	double target = 0.5, constraint = 0.;
	TOTestObjFun2 testObjFun(target, constraint);
	TestRep testTORep;
	testTORep.setRealRep({0.1, 1.});
	std::vector<OutputHandler*> ohVec;
	TopOptChain testOpt(loadOptParamsChain("chaininput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint, target});
	// Test 2nd target
	target = 0.75;
	testObjFun = TOTestObjFun2(target, constraint);
	testOpt = TopOptChain(loadOptParamsChain("chaininput.xml"), &testObjFun, ohVec);
	testTopOpt(testOpt, testTORep, {constraint, target});
	// Test factory creation
	std::unique_ptr<TopOpt> upTO = TopOptFactory::createTopOpt(loadONI("chaininput.xml"), &testObjFun, ohVec);
	testTopOpt(*upTO, testTORep, {constraint, target});
}
