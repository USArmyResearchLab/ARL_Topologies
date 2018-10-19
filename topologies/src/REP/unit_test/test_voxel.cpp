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

#include "voxelrep.h"
#include "torfactory.h"
#include "meshtestns.h"
#include "catch.hpp"
#include <memory>

using namespace Topologies;

// Function for testing various aspects of a VoxelRep
// Defined here so that all construction methods can be tested identically
void testOneVoxelRepInit(TopOptRep& testVox, double reqVolFrac, bool refined = false)
{
	std::size_t rf = refined ? 2 : 1;
	// Check data sizes
	REQUIRE(testVox.getDataSize() == rf*rf*rf*300);
	std::vector<double> realVec;
	testVox.getRealRep(realVec);
	REQUIRE(realVec.size() == rf*rf*rf*300);
	std::vector<int> discVec;
	testVox.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == rf*rf*rf*300);
	std::vector<std::size_t> sizes;
	testVox.getDataSize(sizes);
	REQUIRE(sizes.size() == 3);
	REQUIRE(sizes[0] == rf*10);
	REQUIRE(sizes[1] == rf*2);
	REQUIRE(sizes[2] == rf*15);
	REQUIRE(testVox.getDimension() == 3);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testVox.getDefiningParameters(discreteParams, realParams);
	REQUIRE(discreteParams.size() == 2);
	REQUIRE(discreteParams[0].size() == 5);
	REQUIRE(discreteParams[0][0] == rf*10); // Nx
	REQUIRE(discreteParams[0][1] == rf*2); // Ny
	REQUIRE(discreteParams[0][2] == rf*15); // Nz
	REQUIRE((MeshElementType)discreteParams[0][3] == metHex); // Mesh element type
	REQUIRE((bool)discreteParams[0][4] == false); //useInterpolatoryFilt 
	REQUIRE(discreteParams[1].size() == 6);
	REQUIRE((DimensionType)discreteParams[1][0] == dt3d); // Number of dimensions
	REQUIRE((MeshFileFormat)discreteParams[1][1] == mffStructured); // Mesh file input type
	REQUIRE((UnknownLocation)discreteParams[1][2] == ulElement); // Nodal unknowns
	REQUIRE((FilterType)discreteParams[1][3] == ftNone); // Filter type
	REQUIRE((PenalizationType)discreteParams[1][4] == ptSIMP); // Penalization method
	REQUIRE((ProjectionType)discreteParams[1][5] == ptNone); // Projection method
	REQUIRE(realParams.size() == 3);
	REQUIRE(realParams[0].size() == 5);
	REQUIRE(realParams[0][0] == Approx(0.1)); // Threshold
	REQUIRE(realParams[0][1] == Approx(0.)); // Filter radius
	REQUIRE(realParams[0][2] == Approx(1.)); // width
	REQUIRE(realParams[0][3] == Approx(0.2)); // length
	REQUIRE(realParams[0][4] == Approx(1.5)); // height
	REQUIRE(realParams[1].size() == 2); // Penalization function parameters
	REQUIRE(realParams[1][0] == Approx(2.)); // Penalty power
	REQUIRE(realParams[1][1] == Approx(1e-3)); // Minimum density
	REQUIRE(realParams[2].size() == 0); // Projection function params
	// Check vol frac
	REQUIRE(testVox.computeVolumeFraction() == Approx(reqVolFrac));
	// Check mesh
	std::unique_ptr<TOMesh> chkMesh = testVox.get3DVolumeMesh();
	REQUIRE(chkMesh->getNumNodes() == (rf*10 + 1)*(rf*2 + 1)*(rf*15 + 1));
	REQUIRE(chkMesh->getNumElements() == rf*rf*rf*300);
}

void testOneHeaviRepInit(TopOptRep& testVox, double reqVolFrac, bool refined = false)
{
	std::size_t rf = refined ? 2 : 1;
	// Check data sizes
	REQUIRE(testVox.getDataSize() == (rf*10 + 1)*(rf*2 + 1)*(rf*15 + 1));
	std::vector<double> realVec;
	testVox.getRealRep(realVec);
	REQUIRE(realVec.size() == (rf*10 + 1)*(rf*2 + 1)*(rf*15 + 1));
	std::vector<int> discVec;
	testVox.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == (rf*10 + 1)*(rf*2 + 1)*(rf*15 + 1));
	std::vector<std::size_t> sizes;
	testVox.getDataSize(sizes);
	REQUIRE(sizes.size() == 3);
	REQUIRE(sizes[0] == rf*10 + 1);
	REQUIRE(sizes[1] == rf*2 + 1);
	REQUIRE(sizes[2] == rf*15 + 1);
	REQUIRE(testVox.getDimension() == 3);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testVox.getDefiningParameters(discreteParams, realParams);
	REQUIRE(discreteParams.size() == 2);
	REQUIRE(discreteParams[0].size() == 5);
	REQUIRE(discreteParams[0][0] == rf*10); // Nx
	REQUIRE(discreteParams[0][1] == rf*2); // Ny
	REQUIRE(discreteParams[0][2] == rf*15); // Nz
	REQUIRE((MeshElementType)discreteParams[0][3] == metTet); // Mesh element type
	REQUIRE((bool)discreteParams[0][4] == false); //useInterpolatoryFilt 
	REQUIRE(discreteParams[1].size() == 6);
	REQUIRE((DimensionType)discreteParams[1][0] == dt3d); // Number of dimensions
	REQUIRE((MeshFileFormat)discreteParams[1][1] == mffStructured); // Mesh file input type
	REQUIRE((UnknownLocation)discreteParams[1][2] == ulNode); // Nodal unknowns
	REQUIRE((FilterType)discreteParams[1][3] == ftLinear); // Filter type
	REQUIRE((PenalizationType)discreteParams[1][4] == ptSIMP); // Penalization method
	REQUIRE((ProjectionType)discreteParams[1][5] == ptThresholdHeavi); // Projection method
	REQUIRE(realParams.size() == 3);
	REQUIRE(realParams[0].size() == 5);
	REQUIRE(realParams[0][0] == Approx(0.1)); // Threshold
	REQUIRE(realParams[0][1] == Approx(1.5)); // Filter radius
	REQUIRE(realParams[0][2] == Approx(1.0)); // width
	REQUIRE(realParams[0][3] == Approx(0.2)); // length
	REQUIRE(realParams[0][4] == Approx(1.5)); // height
	REQUIRE(realParams[1].size() == 2); // Penalization function parameters
	REQUIRE(realParams[1][0] == Approx(5.)); // Penalty power
	REQUIRE(realParams[1][1] == Approx(1e-3)); // Minimum density
	REQUIRE(realParams[2].size() == 2); // Projection function params
	REQUIRE(realParams[2][0] == Approx(0.1)); // Threshold
	REQUIRE(realParams[2][1] == Approx(1.e-12)); // Heaviside exponent
	// Check vol frac
	REQUIRE(testVox.computeVolumeFraction() == Approx(reqVolFrac));
	// Check mesh
	std::unique_ptr<TOMesh> chkMesh = testVox.get3DVolumeMesh();
	REQUIRE(chkMesh->getNumNodes() == (rf*10 + 1)*(rf*2 + 1)*(rf*15 + 1));
	REQUIRE(chkMesh->getNumElements() == 6*rf*rf*rf*300); // Triangular elements
}

TEST_CASE("Testing creation of VoxelRep class","[VoxelRep]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testvox.xml");
	std::unique_ptr<TopOptRep> upVox = TopOptRepFactory::createTopOptRep(testRNI);
	testOneVoxelRepInit(*upVox, 1./300.);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upVox->getDefiningParameters(discreteParams, realParams);
	VoxelRep<> testVox(tortVoxel, discreteParams, realParams);
	testVox.initialize(0.234); // Test initialize
	testOneVoxelRepInit(testVox, 0.234);
	// Test assignment operator
	testVox = VoxelRep<>(tortVoxel, discreteParams, realParams);
	std::vector<double> realVec(300, 0.123);
  testVox.setRealRep(realVec.begin(), realVec.end());
	testOneVoxelRepInit(testVox, 0.123);
	// Test copy ctor
	VoxelRep<> copyVox(testVox);
	testOneVoxelRepInit(copyVox, 0.123);
	// Test refine
	copyVox.refine();
	testOneVoxelRepInit(copyVox, 0.123, true);
	// Test MPI rep
	std::vector<std::vector<int>> discreteVars;
	std::vector<std::vector<double>> realVars;
	copyVox.getMPIRep(discreteVars, realVars);
	copyVox.setMPIRep(discreteVars, realVars);
	testOneVoxelRepInit(copyVox, 0.123, true);
	// Test move ctor
	VoxelRep<> moveVox(std::move(copyVox));
	testOneVoxelRepInit(moveVox, 0.123, true);
}

TEST_CASE("Testing creation of Heaviside3D class","[Heaviside3D]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testheavi3.xml");
	std::unique_ptr<TopOptRep> upHeavi = TopOptRepFactory::createTopOptRep(testRNI);
	testOneHeaviRepInit(*upHeavi, 1./300.);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upHeavi->getDefiningParameters(discreteParams, realParams);
	VoxelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> testHeavi(tortHeaviside3D, discreteParams, realParams);
	testHeavi.initialize(0.234); // Test initialize
	testOneHeaviRepInit(testHeavi, 0.234);
	// Test assignment operator
	testHeavi = VoxelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>(tortHeaviside3D, discreteParams, realParams);
	std::vector<double> realVec(528, 0.123);
	testHeavi.setRealRep(realVec.begin(), realVec.end());
	testOneHeaviRepInit(testHeavi, 0.123);
	// Test copy ctor
	VoxelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> copyHeavi(testHeavi);
	testOneHeaviRepInit(copyHeavi, 0.123);
	// Test refine
	copyHeavi.refine();
	testOneHeaviRepInit(copyHeavi, 0.123, true);
	// Test move ctor
	VoxelRep<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> moveHeavi(std::move(copyHeavi));
	testOneHeaviRepInit(moveHeavi, 0.123, true);
}

