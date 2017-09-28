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

#include "genericmaterial.h"
#include "catch.hpp"

using namespace Topologies;

TEST_CASE("Testing default GenericMaterial constructor", "[GenericMaterial]")
{
	GenericMaterial testGM;
	REQUIRE(testGM.getNumParameters() == 0);
	REQUIRE(!testGM.hasLimits());
	std::vector<double> testRGB = {0., 0., 0.};
	std::vector<double> rgb = testGM.getPrintColor();
	REQUIRE(testRGB[0] == rgb[0]);
	REQUIRE(testRGB[1] == rgb[1]);
	REQUIRE(testRGB[2] == rgb[2]);
	GenericMaterial testGM2;
	REQUIRE(testGM.dist(testGM2) == 0.);
	// Copy constructor and operator==
	GenericMaterial testGM3(testGM);
	REQUIRE(testGM == testGM3);
}

TEST_CASE("Testing GenericMaterial constructor with no limits", "[GenericMaterial]")
{
	std::vector<double> matParams = {0., 1., 2.};
	GenericMaterial testGM(matParams);
	REQUIRE(testGM.getNumParameters() == 3);
	for(std::size_t k = 0; k < matParams.size(); ++k)
		REQUIRE(testGM.getParameter(k) == matParams[k]);
	REQUIRE(!testGM.hasLimits());
	std::vector<double> testRGB = {0., 0., 0.};
	std::vector<double> rgb = testGM.getPrintColor();
	REQUIRE(testRGB[0] == rgb[0]);
	REQUIRE(testRGB[1] == rgb[1]);
	REQUIRE(testRGB[2] == rgb[2]);
	GenericMaterial testGM2({0., 0., 0.});
	REQUIRE(testGM.dist(testGM2) == Approx(sqrt(5.)));
	// Copy constructor and operator==
  GenericMaterial testGM3(testGM);
	REQUIRE(testGM == testGM3);
}

TEST_CASE("Testing GenericMaterial constructor with limits", "[GenericMaterial]")
{
	std::vector<double> matParams = {0., -1., 2.};
	std::vector<double> minLim = {0., -1., -2.}, maxLim = {1., 1., 2.};
	GenericMaterial testGM(matParams, minLim, maxLim);
	REQUIRE(testGM.getNumParameters() == 3);
	for(std::size_t k = 0; k < matParams.size(); ++k)
		REQUIRE(testGM.getParameter(k) == matParams[k]);
	REQUIRE(testGM.hasLimits());
	for(std::size_t k = 0; k < matParams.size(); ++k)
	{
		REQUIRE(testGM.getParameterMin(k) == minLim[k]);
		REQUIRE(testGM.getParameterMax(k) == maxLim[k]);
	}
	// print color
	std::vector<double> testRGB = {0., 0., 1.};
	std::vector<double> rgb = testGM.getPrintColor();
	REQUIRE(testRGB[0] == rgb[0]);
	REQUIRE(testRGB[1] == rgb[1]);
	REQUIRE(testRGB[2] == rgb[2]);
	// dist
	GenericMaterial testGM2({0., 0., 0.});
	REQUIRE(testGM.dist(testGM2) == Approx(sqrt(5.)));
	// Property scaling
	testGM.setScaledMaterialParam(0, 0.5);
	REQUIRE(testGM.getParameter(0) == 0.5);
	testGM.setScaledMaterialParam(1, 0.5);
	REQUIRE(testGM.getParameter(1) == 0.);
	testGM.setScaledMaterialParam(2, 0.25);
	REQUIRE(testGM.getParameter(2) == -1.);
	// Random material generation
	std::cout << "random" << std::endl;
	for(unsigned k = 0; k < 100; ++k)
	{
		testGM.genRandomMat();
		// Check that properties are within bounds:
		for(std::size_t kp = 0; kp < matParams.size(); ++kp)
		{
			REQUIRE(testGM.getParameter(kp) <= maxLim[kp]);
			REQUIRE(testGM.getParameter(kp) >= minLim[kp]);
		}
	}
	// Copy constructor and operator==
  GenericMaterial testGM3(testGM);
	REQUIRE(testGM == testGM3);
}

TEST_CASE("Testing GenericMaterial constructor with limits given by a second GenericMaterial", "[GenericMaterial]")
{
	std::vector<double> matParams = {0., -1., 2.};
	std::vector<double> minLim = {0., -1., -2.}, maxLim = {1., 1., 2.};
	GenericMaterial tmpGM(matParams, minLim, maxLim);
	GenericMaterial testGM(matParams, tmpGM);
	REQUIRE(testGM.getNumParameters() == 3);
	for(std::size_t k = 0; k < matParams.size(); ++k)
		REQUIRE(testGM.getParameter(k) == matParams[k]);
	REQUIRE(testGM.hasLimits());
	for(std::size_t k = 0; k < matParams.size(); ++k)
	{
		REQUIRE(testGM.getParameterMin(k) == minLim[k]);
		REQUIRE(testGM.getParameterMax(k) == maxLim[k]);
	}
	// print color
	std::vector<double> testRGB = {0., 0., 1.};
	std::vector<double> rgb = testGM.getPrintColor();
	REQUIRE(testRGB[0] == rgb[0]);
	REQUIRE(testRGB[1] == rgb[1]);
	REQUIRE(testRGB[2] == rgb[2]);
	// dist
	GenericMaterial testGM2({0., 0., 0.});
	REQUIRE(testGM.dist(testGM2) == Approx(sqrt(5.)));
	// Property scaling
	testGM.setScaledMaterialParam(0, 0.5);
	REQUIRE(testGM.getParameter(0) == 0.5);
	testGM.setScaledMaterialParam(1, 0.5);
	REQUIRE(testGM.getParameter(1) == 0.);
	testGM.setScaledMaterialParam(2, 0.25);
	REQUIRE(testGM.getParameter(2) == -1.);
	// Random material generation
	for(unsigned k = 0; k < 100; ++k)
	{
		testGM.genRandomMat();
		// Check that properties are within bounds:
		for(std::size_t kp = 0; kp < matParams.size(); ++kp)
		{
			REQUIRE(testGM.getParameter(kp) <= maxLim[kp]);
			REQUIRE(testGM.getParameter(kp) >= minLim[kp]);
		}
	}
	// Copy constructor and operator==
  GenericMaterial testGM3(testGM);
	REQUIRE(testGM == testGM3);
}

TEST_CASE("Testing GenericMaterial constructor with print color", "[GenericMaterial]")
{
	std::vector<double> matParams = {0., -1., 2.};
  std::vector<double> minLim = {0., -1., -2.}, maxLim = {1., 1., 2.};
  GenericMaterial tmpGM(matParams, minLim, maxLim);
	std::vector<double> rgb = {0., 0.5, 1.};
  GenericMaterial testGM(matParams, tmpGM, rgb);
	std::vector<double> testRGB = testGM.getPrintColor();
	for(std::size_t k = 0; k < matParams.size(); ++k)
		REQUIRE(rgb[k] == testRGB[k]);
	// Copy constructor and operator==
  GenericMaterial testGM3(testGM);
	REQUIRE(testGM == testGM3);
}
