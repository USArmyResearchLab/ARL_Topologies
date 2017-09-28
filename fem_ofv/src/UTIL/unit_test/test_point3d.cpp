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
#include "point3d.h"

TEST_CASE( "Testing Point3D class", "[Point3D]" )
{
	Point3D po(0.,0.,0.);
	Point3D px(1.,0.,0.);
	Point3D py(0.,1.,0.);
	Point3D pz(0.,0.,1.);
	double tol = 1e-14;
	SECTION("Testing point at origin")
	{
		REQUIRE(po.x == 0.);
		REQUIRE(po.y == 0.);
		REQUIRE(po.z == 0.);
		REQUIRE(abs(po) == 0.);
	}
	SECTION("Testing unit vectors")
	{
		REQUIRE(px.x == 1.);
		REQUIRE(px.y == 0.);
		REQUIRE(px.z == 0.);
		REQUIRE(abs(px) == 1.);

		REQUIRE(py.x == 0.);
		REQUIRE(py.y == 1.);
		REQUIRE(py.z == 0.);
		REQUIRE(abs(py) == 1.);

		REQUIRE(pz.x == 0.);
		REQUIRE(pz.y == 0.);
		REQUIRE(pz.z == 1.);
		REQUIRE(abs(pz) == 1.);
	}
	SECTION("Testing logical operators")
	{
		REQUIRE(px == px);
		Point3D t = px;
		REQUIRE(px == t);
	}
	SECTION("Testing arithmetic operators")
	{
		REQUIRE((po + px) == px);

		Point3D t(1.,1.,1.);
		REQUIRE((px + py + pz) == t);

		t = Point3D(2., 2., 2.);
		REQUIRE((2.*px + 2.*py + 2.*pz) == t);
		REQUIRE((px + py + pz) == t/2.);

		po += px;
		po += py;
		po += pz;
		t = Point3D(1.,1.,1.);
		REQUIRE(po == t);

		po = Point3D(0.,0.,0.);
		REQUIRE((-px - py - pz) == -t);
	}
	SECTION("Vector operators")
	{
		// Dot products
		REQUIRE(px*py == 0.);
		REQUIRE(px*pz == 0.);
		REQUIRE(pz*py == 0.);
		REQUIRE(px*px == 1.);
		REQUIRE(py*py == 1.);
		REQUIRE(pz*pz == 1.);
		REQUIRE(px*(2.*py) == 0.);
		REQUIRE(px*(py + 2*px) == 2.);
		REQUIRE(pz*(py + 2*px) == 0.);
		REQUIRE((px + py + pz)*(px + py + pz) == 3.);
		// Cross products
		REQUIRE(px.crossProductNorm(py) == 1.);
		Point3D u = px + py, v = 5.*px + py;
		REQUIRE(u.crossProductNorm(v) == 4.);
		REQUIRE(v.crossProductNorm(u) == 4.);
		REQUIRE(crossProduct(px,py) == pz);
		REQUIRE(crossProduct(py,pz) == px);
		REQUIRE(crossProduct(pz,px) == py);
  }
	SECTION("Area, volume and centroid")
	{
		REQUIRE(triangleArea(po, px, py) == 0.5);
		REQUIRE(triangleArea(po, py, pz) == 0.5);
		REQUIRE(triangleArea(po, pz, px) == 0.5);
		REQUIRE(fabs(tetrahedronVolume(po, px, py, pz) - 1./6.) < tol);
		REQUIRE(abs(centroid(po, px, py, pz) - Point3D(0.25, 0.25, 0.25)) < tol);
	}
	SECTION("Solid angle")
	{
		Point3D p1(1., 0., -1./sqrt(2)), p2(-1., 0., -1./sqrt(2)),
						p3(0., 1., 1./sqrt(2)), p4(0., -1., 1./sqrt(2));
		double exactRes = acos(23./27.);
		REQUIRE(fabs(solidAngle(p1, p2, p3, p4) - exactRes) < tol);
		// Opposite winding
		REQUIRE(fabs(solidAngle(p1, p4, p3, p2) - exactRes) < tol);
	}
}
