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

#include "objfuntest.h"
#include "pixelrep.h"
#include "catch.hpp"
#include <iostream>

using namespace Topologies; 

TEST_CASE("Testing TopOptObjFun finite difference gradient with test derived class", "[TopOptObjFun]")
{
	// Set up PixelRep
	std::vector<std::vector<int>> discreteParams(2);
	int nx = 10, ny = 10;
	discreteParams[0] = {nx, ny, metQuad}; // nx, ny, MeshElementType
	VolMeshTORSpecification tmpVMTORS(tortPixel);
	discreteParams[1] = tmpVMTORS.toVec();
	std::vector<std::vector<double>> realParams(3);
	double width = 1., height = 2.;
	realParams[0] = {0.5, 0., width, height}; // threshold, filterRad, width, height
	realParams[1] = {3., 0.001}; // Penal power, minDensity
	PixelRep<HelperNS::powPenalMin> testPR(tortPixel, discreteParams, realParams);
	testPR.initialize(0.5);
	TOTestObjFun testOF;
	testOF.printResult(testPR, "filename"); // run printResult
	SECTION("f and c")
	{
		std::pair<double, bool> res = testOF(testPR);
		REQUIRE(res.first == Approx(0.5));
		REQUIRE(res.second);
		std::pair<std::vector<double>, bool> outRes;
		testOF.f(testPR, outRes);
		REQUIRE(outRes.first.size() == 1);
		REQUIRE(outRes.first[0] == Approx(0.5));
		REQUIRE(outRes.second);
		testOF.c(testPR, outRes);
		REQUIRE(outRes.first.size() == 1);
		REQUIRE(outRes.first[0] == Approx(0.5));
		REQUIRE(outRes.second);
	}
	SECTION("Gradient of f")
	{
		std::pair<std::vector<double>, bool> outRes;
		testOF.g(testPR, outRes);
		double vol = (1./(double)nx)*(1./(double)ny);
		for(auto oit = outRes.first.begin(); oit != outRes.first.end(); ++oit)
			REQUIRE(*oit == Approx(vol));
	}
	SECTION("Gradient of c")
	{
		std::pair<std::vector<double>, bool> outRes;
		testOF.gc(testPR, outRes);
		double vol = (1./(double)nx)*(1./(double)ny);
		for(auto oit = outRes.first.begin(); oit != outRes.first.end(); ++oit)
			REQUIRE(*oit == Approx(vol));
	}
}
