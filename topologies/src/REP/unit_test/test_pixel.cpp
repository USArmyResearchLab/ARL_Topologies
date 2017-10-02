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

#include "pixelrep.h"
#include "torfactory.h"
#include "meshtestns.h"
#include "catch.hpp"
#include <memory>

using namespace Topologies;

// Function for testing various aspects of a PixelRep
// Defined here so that all construction methods can be tested identically
void testOnePixelRepInit(TopOptRep& testPix, double reqVolFrac, bool refined = false)
{
	std::size_t rf = refined ? 2 : 1;
	// Check data sizes
	REQUIRE(testPix.getDataSize() == rf*rf*450);
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	REQUIRE(realVec.size() == rf*rf*450);
	std::vector<int> discVec;
	testPix.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == rf*rf*450);
	std::vector<std::size_t> sizes;
	testPix.getDataSize(sizes);
	REQUIRE(sizes.size() == 2);
	REQUIRE(sizes[0] == rf*15);
	REQUIRE(sizes[1] == rf*30);
	REQUIRE(testPix.getDimension() == 2);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testPix.getDefiningParameters(discreteParams, realParams);
	REQUIRE(discreteParams.size() == 2);
	REQUIRE(discreteParams[0].size() == 3);
	REQUIRE(discreteParams[0][0] == rf*15); // Nx
	REQUIRE(discreteParams[0][1] == rf*30); // Ny
	REQUIRE((MeshElementType)discreteParams[0][2] == metQuad); // Mesh element type
	REQUIRE(discreteParams[1].size() == 6);
	REQUIRE((DimensionType)discreteParams[1][0] == dt2d); // Number of dimensions
	REQUIRE((MeshFileFormat)discreteParams[1][1] == mffStructured); // Mesh file input type
	REQUIRE((UnknownLocation)discreteParams[1][2] == ulElement); // Nodal unknowns
	REQUIRE((FilterType)discreteParams[1][3] == ftNone); // Filter type
	REQUIRE((PenalizationType)discreteParams[1][4] == ptSIMP); // Penalization method
	REQUIRE((ProjectionType)discreteParams[1][5] == ptNone); // Projection method
	REQUIRE(realParams.size() == 3);
	REQUIRE(realParams[0].size() == 4);
	REQUIRE(realParams[0][0] == Approx(0.5)); // Threshold
	REQUIRE(realParams[0][1] == Approx(0.)); // Filter radius
	REQUIRE(realParams[0][2] == Approx(1.)); // width
	REQUIRE(realParams[0][3] == Approx(2.)); // height
	REQUIRE(realParams[1].size() == 2); // Penalization function parameters
	REQUIRE(realParams[1][0] == Approx(3.)); // Penalty power
	REQUIRE(realParams[1][1] == Approx(1e-3)); // Minimum density
	REQUIRE(realParams[2].size() == 0); // Projection function params
	// Check vol frac
	REQUIRE(testPix.computeVolumeFraction() == Approx(reqVolFrac));
	// Check mesh
	std::unique_ptr<TOMesh> chkMesh = testPix.get2DMesh();
	REQUIRE(chkMesh->getNumNodes() == (rf*15 + 1)*(rf*30 + 1));
	REQUIRE(chkMesh->getNumElements() == rf*rf*450);
}

void testOneHeaviRepInit(TopOptRep& testPix, double reqVolFrac, bool refined = false)
{
	std::size_t rf = refined ? 2 : 1;
	// Check data sizes
	REQUIRE(testPix.getDataSize() == (rf*15 + 1)*(rf*30 + 1));
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	REQUIRE(realVec.size() == (rf*15 + 1)*(rf*30 + 1));
	std::vector<int> discVec;
	testPix.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == (rf*15 + 1)*(rf*30 + 1));
	std::vector<std::size_t> sizes;
	testPix.getDataSize(sizes);
	REQUIRE(sizes.size() == 2);
	REQUIRE(sizes[0] == rf*15 + 1);
	REQUIRE(sizes[1] == rf*30 + 1);
	REQUIRE(testPix.getDimension() == 2);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testPix.getDefiningParameters(discreteParams, realParams);
	REQUIRE(discreteParams.size() == 2);
	REQUIRE(discreteParams[0].size() == 3);
	REQUIRE(discreteParams[0][0] == rf*15); // Nx
	REQUIRE(discreteParams[0][1] == rf*30); // Ny
	REQUIRE((MeshElementType)discreteParams[0][2] == metTri); // Mesh element type
	REQUIRE(discreteParams[1].size() == 6);
  REQUIRE((DimensionType)discreteParams[1][0] == dt2d); // Number of dimensions
  REQUIRE((MeshFileFormat)discreteParams[1][1] == mffStructured); // Mesh file input type
  REQUIRE((UnknownLocation)discreteParams[1][2] == ulNode); // Nodal unknowns
  REQUIRE((FilterType)discreteParams[1][3] == ftLinear); // Filter type
  REQUIRE((PenalizationType)discreteParams[1][4] == ptSIMP); // Penalization method
  REQUIRE((ProjectionType)discreteParams[1][5] == ptThresholdHeavi); // Projection method
	REQUIRE(realParams.size() == 3);
	REQUIRE(realParams[0].size() == 4);
	REQUIRE(realParams[0][0] == Approx(0.5)); // Threshold
	REQUIRE(realParams[0][1] == Approx(1.1)); // Filter radius
	REQUIRE(realParams[0][2] == Approx(1.)); // width
	REQUIRE(realParams[0][3] == Approx(2.)); // height
	REQUIRE(realParams[1].size() == 2); // Penalization function parameters
	REQUIRE(realParams[1][0] == Approx(3.)); // Penalty power
	REQUIRE(realParams[1][1] == Approx(1e-3)); // Minimum density
	REQUIRE(realParams[2].size() == 2); // Projection function params
	REQUIRE(realParams[2][0] == Approx(0.5)); // Threshold
	REQUIRE(realParams[2][1] == Approx(1.e-12)); // Heaviside exponent
	// Check vol frac
	REQUIRE(testPix.computeVolumeFraction() == Approx(reqVolFrac));
	// Check mesh
	std::unique_ptr<TOMesh> chkMesh = testPix.get2DMesh();
	REQUIRE(chkMesh->getNumNodes() == (rf*15 + 1)*(rf*30 + 1));
	REQUIRE(chkMesh->getNumElements() == 2*rf*rf*450); // Triangular elements
}

TEST_CASE("Testing creation of PixelRep class","[PixelRep]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testpix.xml");
	std::unique_ptr<TopOptRep> upPix = TopOptRepFactory::createTopOptRep(testRNI);
	testOnePixelRepInit(*upPix, 1./450.);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upPix->getDefiningParameters(discreteParams, realParams);
	PixelRep<> testPix(tortPixel, discreteParams, realParams);
	testPix.initialize(0.234); // Test initialize
	testOnePixelRepInit(testPix, 0.234);
	// Test assignment operator
	testPix = PixelRep<>(tortPixel, discreteParams, realParams);
	std::vector<double> realVec(450, 0.123);
  testPix.setRealRep(realVec);
	testOnePixelRepInit(testPix, 0.123);
	// Test copy ctor
	PixelRep<> copyPix(testPix);
	testOnePixelRepInit(copyPix, 0.123);
	// Test refine
	copyPix.refine();
	testOnePixelRepInit(copyPix, 0.123, true);
	// Test MPI rep
	std::vector<std::vector<int>> discreteVars;
	std::vector<std::vector<double>> realVars;
	copyPix.getMPIRep(discreteVars, realVars);
	copyPix.setMPIRep(discreteVars, realVars);
	testOnePixelRepInit(copyPix, 0.123, true);
	// Test move ctor
	PixelRep<> movePix(std::move(copyPix));
	testOnePixelRepInit(movePix, 0.123, true);
}

TEST_CASE("Testing creation of Heaviside2D class","[Heaviside2D]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testheavi2.xml");
	std::unique_ptr<TopOptRep> upHeavi = TopOptRepFactory::createTopOptRep(testRNI);
	testOneHeaviRepInit(*upHeavi, 1./450.);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upHeavi->getDefiningParameters(discreteParams, realParams);
	PixelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> testHeavi(tortHeaviside2D, discreteParams, realParams);
	testHeavi.initialize(0.234); // Test initialize
	testOneHeaviRepInit(testHeavi, 0.234);
	// Test assignment operator
	testHeavi = PixelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>(tortHeaviside2D, discreteParams, realParams);
	std::vector<double> realVec(496, 0.123);
	testHeavi.setRealRep(realVec);
	testOneHeaviRepInit(testHeavi, 0.123);
	// Test copy ctor
	PixelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> copyHeavi(testHeavi);
	testOneHeaviRepInit(copyHeavi, 0.123);
	// Test refine
	copyHeavi.refine();
	testOneHeaviRepInit(copyHeavi, 0.123, true);
	// Test move ctor
	PixelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> moveHeavi(std::move(copyHeavi));
	testOneHeaviRepInit(moveHeavi, 0.123, true);
}

