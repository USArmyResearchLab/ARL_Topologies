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

