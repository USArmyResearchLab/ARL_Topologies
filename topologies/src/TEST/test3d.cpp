/*
 *  This file is part of Topologies.
 *
 *  Topologies is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Topologies is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Topologies.  If not, see <http://www.gnu.org/licenses/>.
 */

// Questions?
// Contact: Raymond Wildman, raymond.a.wildman.civ@mail.mil

#define CATCH_CONFIG_MAIN

#include "topoptuniverse.h"
#include "topoptrep.h"
#include "catch.hpp"
#include <vector>
#include <fstream>

using namespace Topologies;

void printRes(const TopOptRep& res, const std::string& fileName)
{
	std::vector<double> realRep;
	res.getRealRep(realRep);
	std::ofstream outFile(fileName.c_str());
	outFile << realRep.size();
	outFile.precision(16);
	for(auto it = realRep.begin(); it != realRep.end(); ++it)
		outFile << " " << *it;
	std::cout << std::endl;
}

std::vector<double> readRes(const std::string& fileName)
{
	std::ifstream inFile(fileName.c_str());
	std::size_t n;
	inFile >> n;
	std::vector<double> outVec(n);
	for(std::size_t k = 0; k < n; ++k)
		inFile >> outVec[k];
	return outVec;
}

TEST_CASE("Testing 3D voxel, element-based densities","[Topologies]")
{
	char progName[] = "test3d";
	char fileName[] = "testvox.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("voxres.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 3D voxel, element-based densities and tetrahedral elements","[Topologies]")
{
	char progName[] = "test3d";
	char fileName[] = "testvoxtet.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("voxrestet.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 3D heaviside with voxel, node-based densities","[Topologies]")
{
	char progName[] = "test2d";
	char fileName[] = "testheavi3.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
//	printRes(*result,"heavi3res.txt"); // Generates result file
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("heavi3res.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 3D tet mesh with a fixed block","[Topologies]")
{
	char progName[] = "test3d";
	char fileName[] = "testmesh3_2blocks.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("mesh3_2blocksres.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 3D tet mesh with a fixed block and Heaviside proj.","[Topologies]")
{
	char progName[] = "test3d";
	char fileName[] = "testheavimesh3_2blocks.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("heavimesh3_2blocksres.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

