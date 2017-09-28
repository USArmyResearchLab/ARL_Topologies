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

#include "volmesh3d.h"
#include "torfactory.h"
#include "meshtestns.h"
#include "catch.hpp"
#include <memory>

using namespace Topologies;

// Function for testing various aspects of a PixelRep
// Defined here so that all construction methods can be tested identically
void testOneVM3DInit(TopOptRep& testPix, double reqVolFrac)
{
	// Check data sizes
	REQUIRE(testPix.getDataSize() == 56);
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	REQUIRE(realVec.size() == 56);
	std::vector<int> discVec;
	testPix.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == 56);
	std::vector<std::size_t> sizes;
	testPix.getDataSize(sizes);
	REQUIRE(sizes.size() == 3);
	REQUIRE(sizes[0] == 56);
	REQUIRE(sizes[1] == 1);
	REQUIRE(sizes[2] == 1);
	REQUIRE(testPix.getDimension() == 3);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testPix.getDefiningParameters(discreteParams, realParams);
	// Discrete
	REQUIRE(discreteParams.size() == 3);
	REQUIRE(discreteParams[0].size() == 12); // File name
	std::string fileNameChk(discreteParams[0].begin(), discreteParams[0].end());
	REQUIRE(fileNameChk == "cylinder.txt");
	REQUIRE(discreteParams[1].size() == 6);
	REQUIRE((DimensionType)discreteParams[1][0] == dt3d); // Number of dimensions
	REQUIRE((MeshFileFormat)discreteParams[1][1] == mffExodus); // Mesh file input type
	REQUIRE((UnknownLocation)discreteParams[1][2] == ulElement); // Nodal unknowns
	REQUIRE((FilterType)discreteParams[1][3] == ftNone); // Filter type
	REQUIRE((PenalizationType)discreteParams[1][4] == ptSIMP); // Penalization method
	REQUIRE((ProjectionType)discreteParams[1][5] == ptNone); // Projection method
	REQUIRE(discreteParams[2].size() == 0); // Fixed blocks
	// Continuous
	REQUIRE(realParams.size() == 4);
	REQUIRE(realParams[0].size() == 8);
	REQUIRE(realParams[0][0] == Approx(0.5)); // Threshold
	REQUIRE(realParams[0][1] == Approx(0.)); // Filter radius, default
	REQUIRE(realParams[0][2] == Approx(1.)); // Tetra size, default
	REQUIRE(realParams[0][3] == Approx(1.)); // Edge size, default
	REQUIRE(realParams[0][4] == Approx(10.)); // Facet angle size, default
	REQUIRE(realParams[0][5] == Approx(1.)); // Facet size size, default
	REQUIRE(realParams[0][6] == Approx(1.)); // Facet distance size, default
	REQUIRE(realParams[0][7] == Approx(1.)); // Cell edge ratio size, default
	REQUIRE(realParams[1].size() == 2); // Penalization function parameters
	REQUIRE(realParams[1][0] == Approx(3.)); // Penalty power
	REQUIRE(realParams[1][1] == Approx(1e-3)); // Minimum density
	REQUIRE(realParams[2].size() == 0); // Projection function parameters
	REQUIRE(realParams[3].size() == 0); // Fixed blocks
	// Check vol frac
	REQUIRE(testPix.computeVolumeFraction() == Approx(reqVolFrac));
	// Check mesh
	std::unique_ptr<TOMesh> chkMesh = testPix.get3DVolumeMesh();
	REQUIRE(chkMesh->getNumNodes() == 111);
	REQUIRE(chkMesh->getNumElements() == 56);
}

void testOneHeaviRepInit(TopOptRep& testPix, double reqVolFrac)
{
	// Check data sizes
	std::unique_ptr<TOMesh> chkMesh = testPix.get3DVolumeMesh(); // Data will match mesh
	REQUIRE(testPix.getDataSize() == chkMesh->getNumNodes());
	std::vector<double> realVec;
	testPix.getRealRep(realVec);
	REQUIRE(realVec.size() == chkMesh->getNumNodes());
	std::vector<int> discVec;
	testPix.getDiscreteRep(discVec);
	REQUIRE(discVec.size() == chkMesh->getNumNodes());
	std::vector<std::size_t> sizes;
	testPix.getDataSize(sizes);
	REQUIRE(sizes.size() == 3);
	REQUIRE(sizes[0] == chkMesh->getNumNodes());
	REQUIRE(sizes[1] == 1);
	REQUIRE(testPix.getDimension() == 3);
	// Input parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	testPix.getDefiningParameters(discreteParams, realParams);
	// Discrete
	REQUIRE(discreteParams.size() == 3);
	REQUIRE(discreteParams[0].size() == 13); // File name
	std::string fileNameChk(discreteParams[0].begin(), discreteParams[0].end());
	REQUIRE(fileNameChk == "rectsolid.stl");
	REQUIRE(discreteParams[1].size() == 6);
	REQUIRE((DimensionType)discreteParams[1][0] == dt3d); // Number of dimensions
	REQUIRE((MeshFileFormat)discreteParams[1][1] == mffSTL); // Mesh file input type
	REQUIRE((UnknownLocation)discreteParams[1][2] == ulNode); // Nodal unknowns
	REQUIRE((FilterType)discreteParams[1][3] == ftLinear); // Filter type
	REQUIRE((PenalizationType)discreteParams[1][4] == ptSIMP); // Penalization method
	REQUIRE((ProjectionType)discreteParams[1][5] == ptThresholdHeavi); // Projection method
	REQUIRE(discreteParams[2].size() == 0); // fixed blocks
	// Continuous
	REQUIRE(realParams.size() == 4);
	REQUIRE(realParams[0].size() == 8);
	REQUIRE(realParams[0][0] == Approx(0.5)); // Threshold
	REQUIRE(realParams[0][1] == Approx(1.1)); // Filter radius
	REQUIRE(realParams[0][2] == Approx(0.5)); // Tetra size
	REQUIRE(realParams[0][3] == Approx(0.5)); // Edge size
	REQUIRE(realParams[0][4] == Approx(25.)); // Facet angle size
	REQUIRE(realParams[0][5] == Approx(0.5)); // Facet size
	REQUIRE(realParams[0][6] == Approx(0.125)); // Facet distance size
	REQUIRE(realParams[0][7] == Approx(3.)); // Cell edge ratio size
	REQUIRE(realParams[1].size() == 2); // Penalization function parameters
	REQUIRE(realParams[1][0] == Approx(3.)); // Penalty power
	REQUIRE(realParams[1][1] == Approx(1e-3)); // Minimum density
	REQUIRE(realParams[2].size() == 2); // Projection parameters
	REQUIRE(realParams[2][0] == Approx(0.5)); // Threshold
	REQUIRE(realParams[2][1] == Approx(1.e-12)); // Heaviside exponent (beta)
	REQUIRE(realParams[3].size() == 0); // Fixed blocks
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

TEST_CASE("Testing creation of VolMesh3D class","[VolMesh3D]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testmesh3d.xml");
	std::unique_ptr<TopOptRep> upVM3D = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upVM3D);
	testOneVM3DInit(*upVM3D, 0.5);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upVM3D->getDefiningParameters(discreteParams, realParams);
	VolMesh3D<> testVM3D(tortMesh3D, discreteParams, realParams);
	testVM3D.initialize(0.234); // Test initialize
	testOneVM3DInit(testVM3D, 0.234);
	// Test assignment operator
	testVM3D = VolMesh3D<>(tortMesh3D, discreteParams, realParams);
	std::vector<double> realVec(56, 0.123);
  testVM3D.setRealRep(realVec);
	testOneVM3DInit(testVM3D, 0.123);
	// Test copy ctor
	VolMesh3D<> copyVM3D(testVM3D);
	testOneVM3DInit(copyVM3D, 0.123);
	// Test MPI rep
	std::vector<std::vector<int>> discreteVars;
	std::vector<std::vector<double>> realVars;
	copyVM3D.getMPIRep(discreteVars, realVars);
	copyVM3D.setMPIRep(discreteVars, realVars);
	testOneVM3DInit(copyVM3D, 0.123);
}

TEST_CASE("Testing creation of HeavisideMesh3D","[HeavisideMesh3D]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testheavimesh3d.xml");
	std::unique_ptr<TopOptRep> upHeavi = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upHeavi);
	std::unique_ptr<TOMesh> upMesh = upHeavi->get3DVolumeMesh();
	INFO("Factory creation");
	testOneHeaviRepInit(*upHeavi, 0.5);
	// Test creation from vectors containing defining parameters
	std::vector<std::vector<int>> discreteParams;
	std::vector<std::vector<double>> realParams;
	upHeavi->getDefiningParameters(discreteParams, realParams);
	VolMesh3D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> testHeavi(tortHeavisideMesh3D, discreteParams, realParams);
	testHeavi.initialize(0.234); // Test initialize
	INFO("Vector ctor");
	testOneHeaviRepInit(testHeavi, 0.234);
	// Test assignment operator
	testHeavi = VolMesh3D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>(tortHeavisideMesh3D, discreteParams, realParams);
	upMesh = testHeavi.get3DVolumeMesh();
	std::vector<double> realVec(upMesh->getNumNodes(), 0.123);
	testHeavi.setRealRep(realVec);
	INFO("Assignment operator and setRealRep");
	testOneHeaviRepInit(testHeavi, 0.123);
	// Test copy ctor
	VolMesh3D<HelperNS::powPenalMin, HelperNS::thresholdHeaviside> copyHeavi(testHeavi);
	INFO("Copy ctor")
	testOneHeaviRepInit(copyHeavi, 0.123);
}

void testDiff(TopOptRep& testTOR)
{
	// Test diffs vs. finite difference approximation
	std::vector<std::map<std::size_t, double>> df = testTOR.diffRep();
	std::vector<double> realVec;
	testTOR.getRealRep(realVec);
	REQUIRE(realVec.size() == df.size());
	std::unique_ptr<TOMesh> upMesh = testTOR.get3DVolumeMesh();
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
		std::unique_ptr<TOMesh> upMeshFD = testTOR.get3DVolumeMesh();
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

TEST_CASE("Testing partial derivatives of VolMesh3D", "[VolMesh3D]")
{
	using namespace MeshTestNS;
	InputLoader::RepNodeInfo testRNI = loadRNI("testmesh3d.xml");
	std::unique_ptr<TopOptRep> upVM3D = TopOptRepFactory::createTopOptRep(testRNI);
	REQUIRE(upVM3D);
	upVM3D->initialize(0.5);
	testDiff(*upVM3D);
	testVFDiff(*upVM3D);
}

TEST_CASE("Testing partial derivatives of Heaviside3D", "[VolMesh3D]")
{
	using namespace MeshTestNS;
	// Test factory creation
	InputLoader::RepNodeInfo testRNI = loadRNI("testheavimesh3d.xml");
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
	InputLoader::RepNodeInfo testRNI = loadRNI("testmesh3d_fixed.xml");
	std::unique_ptr<TopOptRep> upVM2D = TopOptRepFactory::createTopOptRep(testRNI);
	upVM2D->initialize(0.1); // Shouldn't change results
	REQUIRE(upVM2D->computeVolumeFraction() == Approx(1.));
	// Test derivatives
	testDiff(*upVM2D);
	testVFDiff(*upVM2D);
}

