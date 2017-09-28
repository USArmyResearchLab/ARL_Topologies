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

#define CATCH_CONFIG_MAIN

#include "elemedge.h"
#include "point2d.h"
#include <iostream>
#include "catch.hpp"

TEST_CASE( "Testing ElemEdge class", "[ElemEdge]" )
{
	Point2D p0(0.,0.), p1(1.,0.), p2(1., 1.);
	ElemEdge<Point2D> testEdge1(&p0, &p1);
	ElemEdge<Point2D> testEdge2(&p1, &p0);
	ElemEdge<Point2D> testEdge3(&p0, &p2);
	SECTION("Test length")
	{
		double tol = 1e-14;
		REQUIRE(fabs(testEdge1.length() - 1.) < tol);
		REQUIRE(fabs(testEdge2.length() - 1.) < tol);
		REQUIRE(fabs(testEdge3.length() - sqrt(2.)) < tol);
	}
	SECTION("Test hash function")
	{
		Elemedge_hash<Point2D> testhash;
		std::size_t res = testhash(&testEdge1);
		REQUIRE(res == testhash(&testEdge2));
		REQUIRE(res != testhash(&testEdge3));
	}
	SECTION("Test equality")
	{
		Elemedge_equal<Point2D> testequal;
		REQUIRE(testequal(&testEdge1, &testEdge2));
		REQUIRE(!testequal(&testEdge1, &testEdge3));
		REQUIRE(!testequal(&testEdge2, &testEdge3));
		// Alignment
		REQUIRE(testEdge1.isAlignedWith(testEdge2));
		REQUIRE(!testEdge2.isAlignedWith(testEdge2));
		REQUIRE(!testEdge2.isAlignedWith(testEdge3));
	}
	SECTION("Basis functions")
	{
		testEdge1.addNodeBases(testEdge1.node0(), 0);
		testEdge1.addNodeBases(testEdge1.node1(), 1);
		REQUIRE(testEdge1.getGlobalBF(0) == 0);
		REQUIRE(testEdge1.getGlobalBF(1) == 1);
	}
	SECTION("Node check")
	{
		REQUIRE(testEdge1.containsNode(&p0));
		REQUIRE(testEdge1.containsNode(&p1));
		REQUIRE(!testEdge1.containsNode(&p2));
	}
}

