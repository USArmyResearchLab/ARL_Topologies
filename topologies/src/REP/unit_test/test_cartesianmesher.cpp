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

#include "cartesianmesher.h"
#include "catch.hpp"
#include <algorithm>

using namespace Topologies;

TEST_CASE("Testing CartesianMesher in 2D", "[CartesianMesher]")
{
	std::vector<unsigned> ndims = {2, 3};
	std::vector<double> optVals(ndims[0]*ndims[1], 1.);
	do
	{
		// Test tri mesh
		std::unique_ptr<TOMesh> testTriMesh = CartesianMesher::generateMesh(metTri, ndims[0], ndims[1], 1., 1., optVals);
		REQUIRE(testTriMesh);
		REQUIRE(testTriMesh->getNumNodes() == (ndims[0] + 1)*(ndims[1] + 1));
		REQUIRE(testTriMesh->getNumElements() == 2*ndims[0]*ndims[1]);
		// Test quad mesh
		std::unique_ptr<TOMesh> testQuadMesh = CartesianMesher::generateMesh(metQuad, ndims[0], ndims[1], 1., 1., optVals);
		REQUIRE(testQuadMesh);
		REQUIRE(testQuadMesh->getNumNodes() == (ndims[0] + 1)*(ndims[1] + 1));
		REQUIRE(testQuadMesh->getNumElements() == ndims[0]*ndims[1]);
	} while(std::next_permutation(ndims.begin(), ndims.end()));
}

TEST_CASE("Testing CartesianMesher in 3D","[CartesianMesher]")
{
	std::vector<unsigned> ndims = {1, 2, 3};
	std::vector<double> optVals(ndims[0]*ndims[1]*ndims[2], 1.);
	do
	{
		// Test tet mesh
		std::unique_ptr<TOMesh> testTetMesh = CartesianMesher::generateMesh(metTet, ndims[0], ndims[1], ndims[2], 1., 1., 1., optVals);
		REQUIRE(testTetMesh);
		REQUIRE(testTetMesh->getNumNodes() == (ndims[0] + 1)*(ndims[1] + 1)*(ndims[2] + 1));
		REQUIRE(testTetMesh->getNumElements() == 6*ndims[0]*ndims[1]*ndims[2]);
		// Test hex mesh
		std::unique_ptr<TOMesh> testHexMesh = CartesianMesher::generateMesh(metHex, ndims[0], ndims[1], ndims[2], 1., 1., 1., optVals);
		REQUIRE(testHexMesh);
		REQUIRE(testHexMesh->getNumNodes() == (ndims[0] + 1)*(ndims[1] + 1)*(ndims[2] + 1));
		REQUIRE(testHexMesh->getNumElements() == ndims[0]*ndims[1]*ndims[2]);
	} while(std::next_permutation(ndims.begin(), ndims.end()));
}

