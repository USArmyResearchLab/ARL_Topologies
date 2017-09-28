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

#include "csgtree.h"
#include "torfactory.h"
#include "meshtestns.h"
#include "catch.hpp"
#include <memory>

using namespace Topologies;

// Function for testing various aspects of a CSGTreeRep
// Defined here so that all construction methods can be tested identically
void testOneCSGAffine(TopOptRep& testPix, double reqVolFrac)
{
	// Check data sizes
	REQUIRE(testPix.getDataSize() == 3*5);
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	REQUIRE(realVec.size() == 3*5);
	std::vector<int> discVec;
	testPix.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == 0);
	std::vector<std::size_t> sizes;
	testPix.getDataSize(sizes);
	REQUIRE(sizes.size() == 2);
	REQUIRE(sizes[0] == 3*5);
	REQUIRE(sizes[1] == 1);
	REQUIRE(testPix.getDimension() == 2);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testPix.getDefiningParameters(discreteParams, realParams);
	REQUIRE(discreteParams.size() == 1);
	REQUIRE(discreteParams[0].size() == 7);
	REQUIRE(discreteParams[0][0] == 6); // Nx
	REQUIRE(discreteParams[0][1] == 12); // Ny
	REQUIRE(discreteParams[0][2] == 2); // num_shapes x
	REQUIRE(discreteParams[0][3] == 3); // num_shapes y
	REQUIRE(discreteParams[0][4] == 8); // num points per shape
	REQUIRE(discreteParams[0][5] == (int)true); // Use affine
	REQUIRE(discreteParams[0][6] == (int)true); // Shapes are holes
	REQUIRE(realParams.size() == 1);
	REQUIRE(realParams[0].size() == 3);
	REQUIRE(realParams[0][0] == Approx(1.e-3)); // Min density
	REQUIRE(realParams[0][1] == Approx(1.)); // width
	REQUIRE(realParams[0][2] == Approx(2.)); // height
	// Check vol frac
	REQUIRE(testPix.computeVolumeFraction() == Approx(reqVolFrac));
	// Check mesh
	std::unique_ptr<TOMesh> chkMesh = testPix.get2DMesh();
	REQUIRE(chkMesh->getNumNodes() == 259);
	REQUIRE(chkMesh->getNumElements() == 444);
}

void testOneCSGPoint(TopOptRep& testPix, double reqVolFrac)
{
	// Check data sizes
	REQUIRE(testPix.getDataSize() == 2*8*5);
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	REQUIRE(realVec.size() == 2*8*5);
	std::vector<int> discVec;
	testPix.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == 0);
	std::vector<std::size_t> sizes;
	testPix.getDataSize(sizes);
	REQUIRE(sizes.size() == 2);
	REQUIRE(sizes[0] == 2*8*5);
	REQUIRE(sizes[1] == 1);
	REQUIRE(testPix.getDimension() == 2);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testPix.getDefiningParameters(discreteParams, realParams);
	REQUIRE(discreteParams.size() == 1);
	REQUIRE(discreteParams[0].size() == 7);
	REQUIRE(discreteParams[0][0] == 6); // Nx
	REQUIRE(discreteParams[0][1] == 12); // Ny
	REQUIRE(discreteParams[0][2] == 2); // num_shapes x
	REQUIRE(discreteParams[0][3] == 3); // num_shapes y
	REQUIRE(discreteParams[0][4] == 8); // num points per shape
	REQUIRE(discreteParams[0][5] == (int)false); // Use affine
	REQUIRE(discreteParams[0][6] == (int)true); // Shapes are holes
	REQUIRE(realParams.size() == 1);
	REQUIRE(realParams[0].size() == 3);
	REQUIRE(realParams[0][0] == Approx(1.e-3)); // Min density
	REQUIRE(realParams[0][1] == Approx(1.)); // width
	REQUIRE(realParams[0][2] == Approx(2.)); // height
	// Check vol frac
	REQUIRE(testPix.computeVolumeFraction() == Approx(reqVolFrac));
	// Check mesh
	std::unique_ptr<TOMesh> chkMesh = testPix.get2DMesh();
	REQUIRE(chkMesh->getNumNodes() == 259);
	REQUIRE(chkMesh->getNumElements() == 444);
}

TEST_CASE("Testing creation of CSGTreeRep class, using affine opt","[CSGTreeRep]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testcsg.xml");
	std::unique_ptr<TopOptRep> upCSG = TopOptRepFactory::createTopOptRep(testRNI);
	INFO("Factory creation");
	testOneCSGAffine(*upCSG, 0.8833273811);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upCSG->getDefiningParameters(discreteParams, realParams);
	CSGTreeRep testPix(discreteParams, realParams);
	INFO("Creation from vectors");
	testOneCSGAffine(testPix, 0.8833273811);
	// Test assignment operator
	testPix = CSGTreeRep(discreteParams, realParams);
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	realVec[0] = realVec[3] = realVec[6] = realVec[9] = realVec[12] = 1.;
  testPix.setRealRep(realVec);
	INFO("Assignment operator");
	testOneCSGAffine(testPix, 0.2573751256);
	// Test copy ctor
	CSGTreeRep copyPix(testPix);
	INFO("Copy ctor");
	testOneCSGAffine(copyPix, 0.2573751256);
	// Move all polygons
	copyPix.getRealRep(realVec);
	realVec = std::vector<double>(realVec.size(), 1.);
	realVec[0] = realVec[3] = realVec[6] = realVec[9] = realVec[12] = 0.001;
	copyPix.setRealRep(realVec);
	INFO("Moving polygons");
	testOneCSGAffine(copyPix, 1.);
}

TEST_CASE("Testing creation of CSGTreeRep class, using point opt","[CSGTreeRep]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testcsg_point.xml");
	std::unique_ptr<TopOptRep> upCSG = TopOptRepFactory::createTopOptRep(testRNI);
	INFO("Factory creation");
	testOneCSGPoint(*upCSG, 0.8833273811);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upCSG->getDefiningParameters(discreteParams, realParams);
	CSGTreeRep testCSG(discreteParams, realParams);
	INFO("Creation from vectors");
	testOneCSGPoint(testCSG, 0.8833273811);
	// Test assignment operator
	testCSG = CSGTreeRep(discreteParams, realParams);
	std::vector<double> realVec;
  testCSG.getRealRep(realVec);
//	std::cout << "vals = [";
//	for(auto it = realVec.begin(); it != realVec.end(); ++it)
//		std::cout << *it << " ";
//	std::cout << "];" << std::endl;
	// Set nodes to corners of region
	realVec[0] = realVec[1] = realVec[3] = realVec[6] = 0.;
	realVec[2] = realVec[4] = realVec[5] = realVec[7] = 1.;
//	MeshTestNS::printTOMesh(*(testCSG.get2DMesh()));
	testCSG.setRealRep(realVec);
	testCSG.getRealRep(realVec);
//	std::cout << "vals = [";
 // for(auto it = realVec.begin(); it != realVec.end(); ++it)
//    std::cout << *it << " ";
//  std::cout << "];" << std::endl;
//	MeshTestNS::printTOMesh(*(testCSG.get2DMesh()));
	INFO("Assignment operator");
	testOneCSGPoint(testCSG, 0.);
	// Test copy ctor
	CSGTreeRep copyCSG(testCSG);
	INFO("Copy ctor");
	testOneCSGPoint(copyCSG, 0.);
	// Test MPI rep
	std::vector<std::vector<int>> discreteVars;
	std::vector<std::vector<double>> realVars;
	copyCSG.getMPIRep(discreteVars, realVars);
	copyCSG.setMPIRep(discreteVars, realVars);
	INFO("MPI rep");
	testOneCSGPoint(copyCSG, 0.);
}

