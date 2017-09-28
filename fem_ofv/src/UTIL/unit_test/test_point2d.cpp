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
#include "catch.hpp"
#include "point2d.h"

TEST_CASE( "Testing Point2D class", "[Point2D]" )
{
  Point2D po(0.,0.);
  Point2D px(1.,0.);
  Point2D py(0.,1.);
	const double pi = 3.1415926535897932384626433832795;	
  SECTION("Point at origin")
  {
    REQUIRE(po.x == 0.);
    REQUIRE(po.y == 0.);
    REQUIRE(abs(po) == 0.);
  }
  SECTION("Unit vectors")
  {
    REQUIRE(px.x == 1.);
    REQUIRE(px.y == 0.);
    REQUIRE(abs(px) == 1.);
		REQUIRE(arg(px) == 0.);

    REQUIRE(py.x == 0.);
    REQUIRE(py.y == 1.);
    REQUIRE(abs(py) == 1.);
		REQUIRE(fabs(arg(py) - pi/2.) < 1.e-14);
  }
  SECTION("Logical operators")
  {
    REQUIRE(px == px);
    Point2D t = px;
    REQUIRE(px == t);
  }
  SECTION("Arithmetic operators")
  {
    REQUIRE((po + px) == px);

    Point2D t(1.,1.);
    REQUIRE((px + py) == t);

    t = Point2D(2., 2.);
    REQUIRE((2.*px + 2.*py) == t);
    REQUIRE((px + py) == t/2.);

    po += px;
    po += py;
    t = Point2D(1.,1.);
    REQUIRE(po == t);

    po = Point2D(0.,0.);
    REQUIRE((-px - py) == -t);
  }
	SECTION("Vector operators")
	{
		// Dot products
		REQUIRE(px*py == 0.);
		REQUIRE(px*px == 1.);
		REQUIRE(px*(2.*py) == 0.);
		REQUIRE(px*(py + 2*px) == 2.);
		// Cross products
		REQUIRE(px.crossProductNorm(py) == 1.);
		Point2D u = px + py, v = 5.*px + py;
		REQUIRE(u.crossProductNorm(v) == -4.);
		REQUIRE(v.crossProductNorm(u) == 4.);
	}
}
