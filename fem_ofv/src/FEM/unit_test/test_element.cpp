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

#include "lintri.h"
#include "linquad.h"
#include "point2d.h"
#include "UTIL/genericmaterial.h"
#include <iostream>
#include "catch.hpp"

using namespace Topologies;

TEST_CASE( "Testing Element class, with LinearTriangle", "[Element]" )
{
	Point2D p0(0.,0.), p1(1.,0.), p2(0., 1.), p3(1., 1.);
	std::vector<double> mats(3, 1.);
	GenericMaterial gm(mats);
	std::vector<Point2D*> pts1 = {&p0, &p1, &p2};
	LinearTriangle<Point2D> testTri1(pts1, gm);
	std::vector<Point2D*> pts2 = {&p2, &p0, &p1};
	LinearTriangle<Point2D> testTri2(pts2, gm);
	std::vector<Point2D*> pts3 = {&p1, &p2, &p0};
	LinearTriangle<Point2D> testTri3(pts3, gm);
	std::vector<Point2D*> pts4 = {&p0, &p1, &p3};
  LinearTriangle<Point2D> testTri4(pts4, gm);
	SECTION("Test hash function")
	{
		Element_hash<Point2D> testhash;
		std::size_t res = testhash(&testTri1);
		REQUIRE(res == testhash(&testTri2));
		REQUIRE(res == testhash(&testTri3));
		REQUIRE(res != testhash(&testTri4));
	}
	SECTION("Test equality")
	{
		Element_equal<Point2D> testequal;
		REQUIRE(testequal(&testTri1, &testTri2));
		REQUIRE(testequal(&testTri1, &testTri3));
		REQUIRE(testequal(&testTri3, &testTri2));
		REQUIRE(!testequal(&testTri1, &testTri4));
		REQUIRE(!testequal(&testTri2, &testTri4));
		REQUIRE(!testequal(&testTri3, &testTri4));
	}
	SECTION("Basis functions")
	{
		REQUIRE(testTri1.getNumNodes() == 3);
		for(std::size_t k = 0; k < testTri1.getNumNodes(); ++k)
			testTri1.addNodeBases(testTri1.getNode(k), k);
		REQUIRE(testTri1.getGlobalBF(0) == 0);
		REQUIRE(testTri1.getGlobalBF(1) == 1);
		REQUIRE(testTri1.getGlobalBF(2) == 2);
	}
	SECTION("Node check")
	{
		REQUIRE(testTri1.containsNode(&p0));
		REQUIRE(testTri1.containsNode(&p1));
		REQUIRE(testTri1.containsNode(&p2));
		REQUIRE(!testTri1.containsNode(&p3));
	}
}

TEST_CASE( "Testing Element class, with LinQuad", "[Element]" )
{
	Point2D p0(0.,0.), p1(1.,0.), p2(1., 1.), p3(0., 1.), p4(2., 0.);
	std::vector<double> mats(3, 1.);
	GenericMaterial gm(mats);
	std::vector<Point2D*> pts1 = {&p0, &p1, &p2, &p3};
	LinearQuadrilateral<Point2D> testQuad1(pts1, gm);
	std::vector<Point2D*> pts2 = {&p3, &p0, &p1, &p2};
	LinearQuadrilateral<Point2D> testQuad2(pts2, gm);
	std::vector<Point2D*> pts3 = {&p2, &p3, &p0, &p1};
	LinearQuadrilateral<Point2D> testQuad3(pts3, gm);
	std::vector<Point2D*> pts4 = {&p1, &p2, &p3, &p0};
	LinearQuadrilateral<Point2D> testQuad4(pts4, gm);
	std::vector<Point2D*> pts5 = {&p4, &p2, &p3, &p0};
  LinearQuadrilateral<Point2D> testQuad5(pts5, gm);
	SECTION("Test hash struct")
	{
		Element_hash<Point2D> testhash;
		std::size_t res = testhash(&testQuad1);
		REQUIRE(res == testhash(&testQuad2));
		REQUIRE(res == testhash(&testQuad3));
		REQUIRE(res == testhash(&testQuad4));
		REQUIRE(res != testhash(&testQuad5));
	}
	SECTION("Test equality struct")
	{
		Element_equal<Point2D> testequal;
		REQUIRE(testequal(&testQuad1, &testQuad2));
		REQUIRE(testequal(&testQuad1, &testQuad3));
		REQUIRE(testequal(&testQuad1, &testQuad4));
		REQUIRE(testequal(&testQuad2, &testQuad3));
		REQUIRE(testequal(&testQuad2, &testQuad3));
		REQUIRE(testequal(&testQuad3, &testQuad4));
		REQUIRE(!testequal(&testQuad1, &testQuad5));
		REQUIRE(!testequal(&testQuad2, &testQuad5));
		REQUIRE(!testequal(&testQuad3, &testQuad5));
		REQUIRE(!testequal(&testQuad4, &testQuad5));
	}
}
