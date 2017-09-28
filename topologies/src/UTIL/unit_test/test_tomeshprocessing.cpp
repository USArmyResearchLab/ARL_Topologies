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

#include "tomeshprocessing.h"
#include "tomesh.h"
#include "catch.hpp"
#include <algorithm>
#include <iostream>

using namespace Topologies;

TEST_CASE("Testing TOMeshProcessing namespace functions, 2D","[TOMeshProcessing]")
{
	using namespace TOMeshProcessing;
	Point_2_base p0(0.,0.), p1(1.,0.), p2(1.,1.), p3(0.,1.);
	std::vector<Point_2_base> nodeVec = {p0, p1, p2, p3};	
	std::vector<double> optVals(1,1.);
	std::vector<std::vector<std::size_t>> connVec(1);
	SECTION("Quad element, all permutations of node ordering")
	{
		std::vector<std::size_t> elemConn = {0, 1, 2, 3};
		do
		{
			connVec[0] = elemConn;
			TOMesh2D mesh(nodeVec, connVec, optVals);
			REQUIRE(computeElementVolume(0, &mesh) == Approx(1.));
			Point_2_base centroid = getElementCentroid2D(0, &mesh);
			REQUIRE(centroid.x() == Approx(0.5));
			REQUIRE(centroid.y() == Approx(0.5));
		} while(std::next_permutation(elemConn.begin(),elemConn.end()));
	}
	SECTION("Tri element, positive node order")
	{
		std::vector<std::size_t> elemConn = {0, 1, 3};
		do
		{
			connVec[0] = elemConn;
			TOMesh2D mesh(nodeVec, connVec, optVals);
			REQUIRE(computeElementVolume(0, &mesh) == Approx(0.5));
			Point_2_base centroid = getElementCentroid2D(0, &mesh);
			REQUIRE(centroid.x() == Approx(1./3.));
			REQUIRE(centroid.y() == Approx(1./3.));
		} while(std::next_permutation(elemConn.begin(),elemConn.end()));
	}
}

TEST_CASE("Testing TOMeshProcessing namespace functions, 3D","[TOMeshProcessing]")
{
	using namespace TOMeshProcessing;
	Point_3_base p0(0.,0.,0.), p1(1.,0.,0.), p2(1.,1.,0.), p3(0.,1.,0.);
	Point_3_base p4(0.,0.,1.), p5(1.,0.,1.), p6(1.,1.,1.), p7(0.,1.,1.);
	std::vector<Point_3_base> nodeVec = {p0, p1, p2, p3, p4, p5, p6, p7};
	std::vector<double> optVals(1,1.);
	std::vector<std::vector<std::size_t>> connVec(1);
	SECTION("Hex element, all permutations")
	{
		std::vector<std::size_t> elemConn = {0, 1, 2, 3, 4, 5, 6, 7};
		do
		{
			connVec[0] = elemConn;
			TOMesh3D mesh(nodeVec, connVec, optVals);
			REQUIRE(computeElementVolume(0, &mesh) == Approx(1.));
			Point_3_base centroid = getElementCentroid3D(0, &mesh);
			REQUIRE(centroid.x() == Approx(0.5));
			REQUIRE(centroid.y() == Approx(0.5));
			REQUIRE(centroid.z() == Approx(0.5));
		} while(std::next_permutation(elemConn.begin(),elemConn.end()));
	}
	SECTION("Tet element, all permutations")
	{
		std::vector<std::size_t> elemConn = {0, 1, 3, 4};
		do
		{
			connVec[0] = elemConn;
			TOMesh3D mesh(nodeVec, connVec, optVals);
			REQUIRE(computeElementVolume(0, &mesh) == Approx(1./6.));
			Point_3_base centroid = getElementCentroid3D(0, &mesh);
			REQUIRE(centroid.x() == Approx(0.25));
			REQUIRE(centroid.y() == Approx(0.25));
			REQUIRE(centroid.z() == Approx(0.25));
		} while(std::next_permutation(elemConn.begin(),elemConn.end()));
	}
}
