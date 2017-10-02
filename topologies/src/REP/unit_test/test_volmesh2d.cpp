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

#include "volmesh2d.h"
#include "torfactory.h"
#include "meshtestns.h"
#include "catch.hpp"
#include <memory>

using namespace Topologies;

// Function for testing various aspects of a PixelRep
// Defined here so that all construction methods can be tested identically
void testOneVM2DInit(TopOptRep& testPix, double reqVolFrac)
{
	// Check data sizes
	REQUIRE(testPix.getDataSize() == 25);
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	REQUIRE(realVec.size() == 25);
	std::vector<int> discVec;
	testPix.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == 25);
	std::vector<std::size_t> sizes;
	testPix.getDataSize(sizes);
	REQUIRE(sizes.size() == 2);
	REQUIRE(sizes[0] == 25);
	REQUIRE(testPix.getDimension() == 2);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testPix.getDefiningParameters(discreteParams, realParams);
	// Discrete
	REQUIRE(discreteParams.size() == 3);
	REQUIRE(discreteParams[0].size() == 12); // File name
	std::string fileNameChk(discreteParams[0].begin(), discreteParams[0].end());
	REQUIRE(fileNameChk == "testquad.txt");
	REQUIRE(discreteParams[1].size() == 6);
	REQUIRE((DimensionType)discreteParams[1][0] == dt2d); // Number of dimensions
	REQUIRE((MeshFileFormat)discreteParams[1][1] == mffExodus); // Mesh file input type
	REQUIRE((UnknownLocation)discreteParams[1][2] == ulElement); // Nodal unknowns
	REQUIRE((FilterType)discreteParams[1][3] == ftNone); // Filter type
	REQUIRE((PenalizationType)discreteParams[1][4] == ptSIMP); // Penalization method
	REQUIRE((ProjectionType)discreteParams[1][5] == ptNone); // Projection method
	REQUIRE(discreteParams[2].size() == 0); // Fixed block info
	// Continuous
	REQUIRE(realParams.size() == 5);
	REQUIRE(realParams[0].size() == 4);
	REQUIRE(realParams[0][0] == Approx(0.25)); // Threshold
	REQUIRE(realParams[0][1] == Approx(0.)); // Filter radius, default
	REQUIRE(realParams[0][2] == Approx(10.)); // Mesh element angle, default
	REQUIRE(realParams[0][3] == Approx(1.)); // Mesh edge size, default
	REQUIRE(realParams[1].size() == 2); // Penalization function parameters
	REQUIRE(realParams[1][0] == Approx(4.)); // Penalty power
	REQUIRE(realParams[1][1] == Approx(1e-3)); // Minimum density
	REQUIRE(realParams[2].size() == 0); // Projection function parameters
	REQUIRE(realParams[3].size() == 0); // Fixed block info
	REQUIRE(realParams[4].size() == 40); // Boundary edges
	// Check vol frac
	REQUIRE(testPix.computeVolumeFraction() == Approx(reqVolFrac));
	// Check mesh
	std::unique_ptr<TOMesh> chkMesh = testPix.get2DMesh();
	REQUIRE(chkMesh->getNumNodes() == 36);
	REQUIRE(chkMesh->getNumElements() == 25);
}

void testOneHeaviRepInit(TopOptRep& testPix, double reqVolFrac)
{
	// Check data sizes
	std::unique_ptr<TOMesh> chkMesh = testPix.get2DMesh(); // Data will match mesh
	REQUIRE(testPix.getDataSize() == chkMesh->getNumNodes());
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	REQUIRE(realVec.size() == chkMesh->getNumNodes());
	std::vector<int> discVec;
	testPix.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == chkMesh->getNumNodes());
	std::vector<std::size_t> sizes;
	testPix.getDataSize(sizes);
	REQUIRE(sizes.size() == 2);
	REQUIRE(sizes[0] == chkMesh->getNumNodes());
	REQUIRE(testPix.getDimension() == 2);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testPix.getDefiningParameters(discreteParams, realParams);
	// Discrete
	REQUIRE(discreteParams.size() == 3);
	REQUIRE(discreteParams[0].size() == 2); // File name
	std::string fileNameChk(discreteParams[0].begin(), discreteParams[0].end());
	REQUIRE(fileNameChk == "na");
	REQUIRE(discreteParams[1].size() == 6);
	REQUIRE((DimensionType)discreteParams[1][0] == dt2d); // Number of dimensions
	REQUIRE((MeshFileFormat)discreteParams[1][1] == mffPolygon); // Mesh file input type
	REQUIRE((UnknownLocation)discreteParams[1][2] == ulNode); // Nodal unknowns
	REQUIRE((FilterType)discreteParams[1][3] == ftLinear); // Filter type
	REQUIRE((PenalizationType)discreteParams[1][4] == ptSIMP); // Penalization method
	REQUIRE((ProjectionType)discreteParams[1][5] == ptThresholdHeavi); // Projection method
	REQUIRE(discreteParams[2].size() == 0);
	// Continuous
	REQUIRE(realParams.size() == 6);
	REQUIRE(realParams[0].size() == 4);
	REQUIRE(realParams[0][0] == Approx(0.75)); // Threshold
	REQUIRE(realParams[0][1] == Approx(2.2)); // Filter radius, default
	REQUIRE(realParams[0][2] == Approx(15.)); // Mesh element angle
	REQUIRE(realParams[0][3] == Approx(0.1)); // Mesh edge size
	REQUIRE(realParams[1].size() == 2); // Penalization function parameters
	REQUIRE(realParams[1][0] == Approx(2.)); // Penalty power
	REQUIRE(realParams[1][1] == Approx(1e-3)); // Minimum density
	REQUIRE(realParams[2].size() == 2); // Projection parameters
	REQUIRE(realParams[2][0] == Approx(0.75)); // Threshold
	REQUIRE(realParams[2][1] == Approx(1.e-12)); // Heaviside exponent (beta)
	REQUIRE(realParams[3].size() == 0); // Fixed block info
	REQUIRE(realParams[4].size() == 6); // Boundary edges
	REQUIRE(realParams[5].size() == 8); // Boundary edges
	// Check vol frac
	REQUIRE(testPix.computeVolumeFraction() == Approx(reqVolFrac));
}

double computeEps()
{
	double mEps = 1.;
	do
		mEps /= 2.;
	while (1. + mEps > 1.);
	mEps *= 2.;
	return mEps;
}

TEST_CASE("Testing creation of VolMesh2D class","[VolMesh2D]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testmesh2d.xml");
	std::unique_ptr<TopOptRep> upVM2D = TopOptRepFactory::createTopOptRep(testRNI);
	testOneVM2DInit(*upVM2D, 0.25);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upVM2D->getDefiningParameters(discreteParams, realParams);
	VolMesh2D<> testVM2D(tortMesh2D, discreteParams, realParams);
	testVM2D.initialize(0.234); // Test initialize
	testOneVM2DInit(testVM2D, 0.234);
	// Test assignment operator
	testVM2D = VolMesh2D<>(tortMesh2D, discreteParams, realParams);
	std::vector<double> realVec(25, 0.123);
  testVM2D.setRealRep(realVec);
	testOneVM2DInit(testVM2D, 0.123);
	// Test copy ctor
	VolMesh2D<> copyVM2D(testVM2D);
	testOneVM2DInit(copyVM2D, 0.123);
	// Test MPI rep
	std::vector<std::vector<int>> discreteVars;
	std::vector<std::vector<double>> realVars;
	copyVM2D.getMPIRep(discreteVars, realVars);
	copyVM2D.setMPIRep(discreteVars, realVars);
	testOneVM2DInit(copyVM2D, 0.123);
	// Test move ctor
	VolMesh2D<> moveVM2D(std::move(copyVM2D));
	testOneVM2DInit(moveVM2D, 0.123);
}

TEST_CASE("Testing creation of HeavisideMesh2D","[Heaviside2D]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testheavimesh2d.xml");
	std::unique_ptr<TopOptRep> upHeavi = TopOptRepFactory::createTopOptRep(testRNI);
	INFO("Factory creation");
	testOneHeaviRepInit(*upHeavi, 0.75);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upHeavi->getDefiningParameters(discreteParams, realParams);
	VolMesh2D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> testHeavi(tortHeavisideMesh2D, discreteParams, realParams);
	testHeavi.initialize(0.234); // Test initialize
	INFO("Creation from vectors");
	testOneHeaviRepInit(testHeavi, 0.234);
	// Test assignment operator
	testHeavi = VolMesh2D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>(tortHeavisideMesh2D, discreteParams, realParams);
	std::unique_ptr<TOMesh> upMesh = testHeavi.get2DMesh();
	std::vector<double> realVec(upMesh->getNumNodes(), 0.123);
	testHeavi.setRealRep(realVec);
	INFO("Assignment operator");
	testOneHeaviRepInit(testHeavi, 0.123);
	// Test copy ctor
	VolMesh2D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> copyHeavi(testHeavi);
	INFO("Copy ctor");
	testOneHeaviRepInit(copyHeavi, 0.123);
	// Test move ctor
	VolMesh2D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> moveHeavi(std::move(copyHeavi));
	INFO("Move ctor");
	testOneHeaviRepInit(moveHeavi, 0.123);
}

void testDiff(TopOptRep& testTOR)
{
	// Test diffs vs. finite difference approximation
	std::vector<std::map<std::size_t, double>> df = testTOR.diffRep();
	std::vector<double> realVec;
	testTOR.getRealRep(realVec);
	REQUIRE(realVec.size() == df.size());
	std::unique_ptr<TOMesh> upMesh = testTOR.get2DMesh();
	REQUIRE(upMesh);
	double sqrtEta = sqrt(computeEps());
	for(std::size_t k = 0; k < realVec.size(); ++k)
	{
		// Perturb kth value of realVec
		double tmpk = realVec[k];
		realVec[k] += sqrtEta;
		double step = realVec[k] - tmpk;
		testTOR.setRealRep(realVec);
		// Get new mesh
		std::unique_ptr<TOMesh> upMeshFD = testTOR.get2DMesh();
		REQUIRE(upMeshFD);
		// Each mesh element should match the values in df:
		for(std::size_t ke = 0; ke < upMesh->getNumElements(); ++ke)
		{
			auto it = df[k].find(ke);
			double exact = it != df[k].end() ? it->second : 0.;
			double fdapprox = (upMeshFD->getOptVal(ke) - upMesh->getOptVal(ke))/step;
			REQUIRE(exact == Approx(fdapprox));
		}
		// Set back to original value
		realVec[k] = tmpk;
		testTOR.setRealRep(realVec);
	}
}

void testVFDiff(TopOptRep& testTOR)
{
	// Test vf diffs vs. finite difference approximation
	std::vector<double> vfGrad = testTOR.computeGradVolumeFraction();
	std::vector<double> realVec;
	testTOR.getRealRep(realVec);
	REQUIRE(realVec.size() == vfGrad.size());
	double vf0 = testTOR.computeVolumeFraction();
	double sqrtEta = sqrt(computeEps());
	for(std::size_t k = 0; k < realVec.size(); ++k)
	{
		// Perturb kth value of realVec
		double tmpk = realVec[k];
		realVec[k] += sqrtEta;
		double step = realVec[k] - tmpk;
		testTOR.setRealRep(realVec);
		double vf1 = testTOR.computeVolumeFraction();
		// Compute estimate and compare
		double fdapprox = (vf1 - vf0)/step;
		REQUIRE(vfGrad[k] == Approx(fdapprox));
		// Set back to original value
		realVec[k] = tmpk;
		testTOR.setRealRep(realVec);
	}
}

TEST_CASE("Testing partial derivatives of VolMesh2D", "[VolMesh2D]")
{
	using namespace MeshTestNS;
	InputLoader::RepNodeInfo testRNI = loadRNI("testmesh2d.xml");
	std::unique_ptr<TopOptRep> upVM2D = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upVM2D);
	upVM2D->initialize(0.5);
	testDiff(*upVM2D);
	testVFDiff(*upVM2D);
}

TEST_CASE("Testing partial derivatives of Heaviside2D", "[VolMesh2D]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testheavimesh2d.xml");
	std::unique_ptr<TopOptRep> upHeavi = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upHeavi);
	upHeavi->initialize(0.5);
	testDiff(*upHeavi);
	testVFDiff(*upHeavi);
}

TEST_CASE("Testing VolMesh2D with a fixed block","[VolMesh2D]")
{
  using namespace MeshTestNS;
  // Test factory creation
  InputLoader::RepNodeInfo testRNI = loadRNI("testmesh2d_fixed.xml");
  std::unique_ptr<TopOptRep> upVM2D = TopOptRepFactory::createTopOptRep(testRNI);
	upVM2D->initialize(0.1); // Shouldn't change results
	REQUIRE(upVM2D->computeVolumeFraction() == Approx(1.));
	// Test derivatives
	testDiff(*upVM2D);
  testVFDiff(*upVM2D);
}

