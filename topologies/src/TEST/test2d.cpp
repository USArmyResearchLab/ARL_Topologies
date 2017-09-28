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

TEST_CASE("Testing 2D pixel, element-based densities","[Topologies]")
{
	char progName[] = "test2d";
	char fileName[] = "testpix.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("pixres.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 2D pixel, element-based densities with triangular mesh","[Topologies]")
{
	char progName[] = "test2d";
	char fileName[] = "testpixtri.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("pixrestri.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 2D heaviside with pixel, node-based densities","[Topologies]")
{
	char progName[] = "test2d";
	char fileName[] = "testheavi2.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("heavi2res.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 2D mesh, auto-generated mesh","[Topologies]")
{
	char progName[] = "test2d";
	char fileName[] = "testmesh2auto.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("mesh2autores.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 2D mesh, exodus mesh with fixed block","[Topologies]")
{
	char progName[] = "test2d";
	char fileName[] = "testmesh_2blocks.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("mesh2blockres.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 2D heaviside with exodus mesh and fixed block","[Topologies]")
{
	char progName[] = "test2d";
	char fileName[] = "testheavimesh2_2blocks.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("heavimesh2blockres.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

TEST_CASE("Testing 2D mesh, gmsh mesh with fixed block","[Topologies]")
{
	char progName[] = "test2d";
	char fileName[] = "testgmsh_2blocks.xml";
	char* argv[] = {progName, fileName};
	TopOptUniverse toProblem(2, argv);
	toProblem.runProblem();
	std::unique_ptr<TopOptRep> result(toProblem.getResult());
//	printRes(*result, "gmsh2blockres.txt");
	// Load result from file (generated previously)
	std::vector<double> resVec = readRes("gmsh2blockres.txt");
	// Get optimization values
	std::vector<double> realRep;
	result->getRealRep(realRep);
	// Compare
	REQUIRE(realRep.size() == resVec.size());
	for(std::size_t k = 0; k < realRep.size(); ++k)
		REQUIRE(realRep[k] == Approx(resVec[k]));
}

