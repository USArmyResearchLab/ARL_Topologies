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
#include "coordinatesystem.h"
#include "point2d.h"
#include "point3d.h"

TEST_CASE( "Testing cylindrical coordinates in 2d", "[CoordinateSystem]" )
{
	using namespace CoordinateSystem;
	double const pi = 3.14159265358979323846;
	double const dth = 2.*pi/8.;
	unsigned kth = 0;
	std::vector<Point2D> vVec(8);
	for(auto& v : vVec)
		v = Point2D(1., static_cast<double>(kth++)*dth);
	SECTION("Point conversion")
	{
		// Around a unit circle (constant r), 45 degrees at a time
		kth = 0;
		// 0 degrees
		Point2D res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(1.));
		REQUIRE(res.y == Approx(0.));
		// 45 degrees
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(sqrt(2.)*0.5));
		REQUIRE(res.y == Approx(sqrt(2.)*0.5));
		// 90
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(1.));
		// 135
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(-sqrt(2.)*0.5));
		REQUIRE(res.y == Approx(sqrt(2.)*0.5));
		// 180
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(-1.));
		REQUIRE(res.y == Approx(0.));
		// 225
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(-sqrt(2.)*0.5));
		REQUIRE(res.y == Approx(-sqrt(2.)*0.5));
		// 270
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(-1.));
		// 315
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(sqrt(2.)*0.5));
		REQUIRE(res.y == Approx(-sqrt(2.)*0.5));
		// Along a ray (constant theta)
		res = cylToCart(Point2D(1., pi/3.));
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(sqrt(3.)*0.5));
		res = cylToCart(Point2D(2., pi/3.));
		REQUIRE(res.x == Approx(1.));
		REQUIRE(res.y == Approx(sqrt(3.)));
		res = cylToCart(Point2D(10., pi/3.));
		REQUIRE(res.x == Approx(5.));
		REQUIRE(res.y == Approx(sqrt(3.)*5.));
	}
	SECTION("Vector conversion")
	{
		// Test vectors at various points on the unit circle
		kth = 0;
		for(kth = 0; kth < 8; ++kth)
		{
			Point2D p1 = cylToCart(vVec[kth]);
			Point2D p2 = cylToCart(vVec[kth]); // Results should be the same
			p2 *= pi;
			// Cycle through vectors
			for(unsigned kv = 0; kv < 8; ++kv)
			{
				// p1
				Point2D v = cylToCart(vVec[kv]); // Use existing points as vectors
				Point2D res = cylVectorToCartVector(v, p1);
				Point2D ans = cylToCart(vVec[(kv + kth)%8]);
				REQUIRE(res.x == Approx(ans.x));
				REQUIRE(res.y == Approx(ans.y));
				res = cylVectorToCartVector(v, p2);
				REQUIRE(res.x == Approx(ans.x));
				REQUIRE(res.y == Approx(ans.y));
			}
		}
	}
}

TEST_CASE( "Testing cylindrical coordinates in 3d", "[CoordinateSystem]" )
{
	using namespace CoordinateSystem;
	double const pi = 3.14159265358979323846;
	double const dth = 2.*pi/8.;
	unsigned kth = 0;
	std::vector<Point3D> vVec(8);
	for(auto& v : vVec)
		v = Point3D(1., static_cast<double>(kth++)*dth, 2.);
	SECTION("Point conversion")
	{
		// Around a unit circle (constant r and z), 45 degrees at a time
		kth = 0;
		// 0 degrees
		Point3D res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(1.));
		REQUIRE(res.y == Approx(0.));
		REQUIRE(res.z == Approx(2.));
		// 45 degrees
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(sqrt(2.)*0.5));
		REQUIRE(res.y == Approx(sqrt(2.)*0.5));
		REQUIRE(res.z == Approx(2.));
		// 90
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(1.));
		REQUIRE(res.z == Approx(2.));
		// 135
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(-sqrt(2.)*0.5));
		REQUIRE(res.y == Approx(sqrt(2.)*0.5));
		REQUIRE(res.z == Approx(2.));
		// 180
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(-1.));
		REQUIRE(res.y == Approx(0.));
		REQUIRE(res.z == Approx(2.));
		// 225
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(-sqrt(2.)*0.5));
		REQUIRE(res.y == Approx(-sqrt(2.)*0.5));
		REQUIRE(res.z == Approx(2.));
		// 270
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(-1.));
		REQUIRE(res.z == Approx(2.));
		// 315
		res = cylToCart(vVec[kth++]);
		REQUIRE(res.x == Approx(sqrt(2.)*0.5));
		REQUIRE(res.y == Approx(-sqrt(2.)*0.5));
		REQUIRE(res.z == Approx(2.));
		// Along a ray (constant theta and z)
		res = cylToCart(Point3D(1., pi/3., 2.));
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(sqrt(3.)*0.5));
		REQUIRE(res.z == Approx(2.));
		res = cylToCart(Point3D(2., pi/3., 2.));
		REQUIRE(res.x == Approx(1.));
		REQUIRE(res.y == Approx(sqrt(3.)));
		REQUIRE(res.z == Approx(2.));
		res = cylToCart(Point3D(10., pi/3., 2.));
		REQUIRE(res.x == Approx(5.));
		REQUIRE(res.y == Approx(sqrt(3.)*5.));
		REQUIRE(res.z == Approx(2.));
		// Along height (constant r and theta)
		res = cylToCart(Point3D(1., pi/3., 0.5));
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(sqrt(3.)*0.5));
		REQUIRE(res.z == Approx(0.5));
		res = cylToCart(Point3D(1., pi/3., 5.));
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(sqrt(3.)*0.5));
		REQUIRE(res.z == Approx(5.));
	}
	SECTION("Vector conversion")
	{
		// Test vectors at various points on the unit circle
		kth = 0;
		for(kth = 0; kth < 8; ++kth)
		{
			Point3D p1 = cylToCart(vVec[kth]);
			Point3D p2 = cylToCart(vVec[kth]); // Results should be the same
			p2 *= pi;
			// Cycle through vectors
			for(unsigned kv = 0; kv < 8; ++kv)
			{
				// p1
				Point3D v = cylToCart(vVec[kv]); // Use existing points as vectors
				Point3D res = cylVectorToCartVector(v, p1);
				Point3D ans = cylToCart(vVec[(kv + kth)%8]);
				REQUIRE(res.x == Approx(ans.x));
				REQUIRE(res.y == Approx(ans.y));
				REQUIRE(res.z == Approx(2.));
				res = cylVectorToCartVector(v, p2);
				REQUIRE(res.x == Approx(ans.x));
				REQUIRE(res.y == Approx(ans.y));
				REQUIRE(res.z == Approx(2.));
				// Test generic function
				res = convertVector(v, p2, Type::cylindrical);
				REQUIRE(res.x == Approx(ans.x));
				REQUIRE(res.y == Approx(ans.y));
				REQUIRE(res.z == Approx(2.));
			}
		}
	}
}

TEST_CASE( "Testing spherical coordinates", "[CoordinateSystem]" )
{
	using namespace CoordinateSystem;
	double const pi = 3.14159265358979323846;
	SECTION("Point conversion")
	{
		// Test a point in each octant
		Point3D res = sphrToCart(Point3D(1., pi*0.25, pi*0.25));
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(0.5));
		REQUIRE(res.z == Approx(sqrt(2.)/2.));
		res = sphrToCart(Point3D(1., pi*0.25, pi*0.75));
		REQUIRE(res.x == Approx(-0.5));
		REQUIRE(res.y == Approx(0.5));
		REQUIRE(res.z == Approx(sqrt(2.)/2.));
		res = sphrToCart(Point3D(1., pi*0.25, pi*1.25));
		REQUIRE(res.x == Approx(-0.5));
		REQUIRE(res.y == Approx(-0.5));
		REQUIRE(res.z == Approx(sqrt(2.)/2.));
		res = sphrToCart(Point3D(1., pi*0.25, pi*1.75));
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(-0.5));
		REQUIRE(res.z == Approx(sqrt(2.)/2.));
		res = sphrToCart(Point3D(1., pi*0.75, pi*0.25));
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(0.5));
		REQUIRE(res.z == Approx(-sqrt(2.)/2.));
		res = sphrToCart(Point3D(1., pi*0.75, pi*0.75));
		REQUIRE(res.x == Approx(-0.5));
		REQUIRE(res.y == Approx(0.5));
		REQUIRE(res.z == Approx(-sqrt(2.)/2.));
		res = sphrToCart(Point3D(1., pi*0.75, pi*1.25));
		REQUIRE(res.x == Approx(-0.5));
		REQUIRE(res.y == Approx(-0.5));
		REQUIRE(res.z == Approx(-sqrt(2.)/2.));
		res = sphrToCart(Point3D(1., pi*0.75, pi*1.75));
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(-0.5));
		REQUIRE(res.z == Approx(-sqrt(2.)/2.));
	}
	SECTION("Vector conversion")
	{
		// Test x axis
		Point3D p = sphrToCart(Point3D(1., pi*0.5, 0.));
		Point3D v(1., 0., 0.); // +r
		Point3D res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(p.x));
		REQUIRE(res.y == Approx(p.y));
		REQUIRE(res.z == Approx(p.z));
		v = Point3D(-1., 0., 0.); // -r
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(-p.x));
		REQUIRE(res.y == Approx(-p.y));
		REQUIRE(res.z == Approx(-p.z));
		v = Point3D(0., 1., 0.); // +theta
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(0.));
		REQUIRE(res.z == Approx(-1.));
		v = Point3D(0., -1., 0.); // -theta
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(0.));
		REQUIRE(res.z == Approx(1.));
		v = Point3D(0., 0., 1.); // +phi
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(1.));
		REQUIRE(res.z == Approx(0.));
		v = Point3D(0., 0., -1.); // -phi
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(-1.));
		REQUIRE(res.z == Approx(0.));
		// Test +++ octant
		p = sphrToCart(Point3D(1., pi*0.25, pi*0.25));
		v = Point3D(1., 0., 0.);
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(p.x));
		REQUIRE(res.y == Approx(p.y));
		REQUIRE(res.z == Approx(p.z));
		v = Point3D(-1., 0., 0.);
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(-p.x));
		REQUIRE(res.y == Approx(-p.y));
		REQUIRE(res.z == Approx(-p.z));
		v = Point3D(0., 1., 0.);
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(0.5));
		REQUIRE(res.y == Approx(0.5));
		REQUIRE(res.z == Approx(-sqrt(2.)/2.));
		v = Point3D(0., -1., 0.);
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(-0.5));
		REQUIRE(res.y == Approx(-0.5));
		REQUIRE(res.z == Approx(sqrt(2.)/2.));
		v = Point3D(0., 0., 1.);
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(-sqrt(2.)/2.));
		REQUIRE(res.y == Approx(sqrt(2.)/2.));
		REQUIRE(res.z == Approx(0.));
		v = Point3D(0., 0., -1.);
		res = sphrVectorToCartVector(v, p);
		REQUIRE(res.x == Approx(sqrt(2.)/2.));
		REQUIRE(res.y == Approx(-sqrt(2.)/2.));
		REQUIRE(res.z == Approx(0.));
		// Test generic function on z-axis
		p = sphrToCart(Point3D(2., 0., 0.));
		v = Point3D(1., 0., 0.);
		res = convertVector(v, p, Type::spherical);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(0.));
		REQUIRE(res.z == Approx(1.));
		v = Point3D(-1., 0., 0.);
		res = convertVector(v, p, Type::spherical);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(0.));
		REQUIRE(res.z == Approx(-1.));
		v = Point3D(0., 1., 0.);
		res = convertVector(v, p, Type::spherical);
		REQUIRE(res.x == Approx(1.));
		REQUIRE(res.y == Approx(0.));
		REQUIRE(res.z == Approx(0.));
		v = Point3D(0., -1., 0.);
		res = convertVector(v, p, Type::spherical);
		REQUIRE(res.x == Approx(-1.));
		REQUIRE(res.y == Approx(0.));
		REQUIRE(res.z == Approx(0.));
		v = Point3D(0., 0., 1.);
		res = convertVector(v, p, Type::spherical);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(1.));
		REQUIRE(res.z == Approx(0.));
		v = Point3D(0., 0., -1.);
		res = convertVector(v, p, Type::spherical);
		REQUIRE(res.x == Approx(0.));
		REQUIRE(res.y == Approx(-1.));
		REQUIRE(res.z == Approx(0.));

	}
}

