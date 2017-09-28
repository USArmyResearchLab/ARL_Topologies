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

#include "geometricentity.h"

using namespace Topologies;

TEST_CASE( "Testing GeometricEntity::Point class", "[GeometricEntity::Point]" )
{
	Point_2_base po2(0., 0.);
	Point_3_base po3(0., 0., 0.);
	Point geP2(po2), geP3(po3);
	double tol = 1e-14;
	SECTION("2D point")
	{
		// Point coincident
		REQUIRE(geP2.isPointCoincident(po2, tol));
		// Line segement
		Point_2_base p0(-1, 0.), p1(1.,0.), p2(0., -1.), p3(0.,1.);
		REQUIRE(geP2.doesLineSegmentIntersect(p0, p1));
		REQUIRE(geP2.doesLineSegmentIntersect(p2, p3));
		REQUIRE(!geP2.doesLineSegmentIntersect(p0, p2));
		REQUIRE(!geP2.doesLineSegmentIntersect(p0, p3));
	}
	SECTION("3D point")
	{
		// Point coincident
		REQUIRE(geP3.isPointCoincident(po3, tol));
		// Line segement
		Point_3_base p0(-1, 0.,0.), p1(1.,0.,0.), p2(0.,-1.,0.), p3(0.,1.,0.), p4(0.,0.,-1.), p5(0.,0.,1.);
		REQUIRE(geP3.doesLineSegmentIntersect(p0, p1));
		REQUIRE(geP3.doesLineSegmentIntersect(p2, p3));
		REQUIRE(geP3.doesLineSegmentIntersect(p4, p5));
		REQUIRE(!geP3.doesLineSegmentIntersect(p0, p2));
		REQUIRE(!geP3.doesLineSegmentIntersect(p0, p3));
		REQUIRE(!geP3.doesLineSegmentIntersect(p0, p5));
	}
}

TEST_CASE("Testing GeometricEntity::InfiniteLine class", "[GeometricEntity::InfiniteLine]")
{
	const double pi = 3.1415926535897932384626433832795;
	Point_2_base po(0., 0.), px(1.,0.), py(0.,1.), mpx(-1.,0.);
	double angle = 0.;
	double m = 0., b = 0.;
	double intercept = 1.;
	double tol = 1e-14;
	Point_2_base u(1., 1.), v(-1., -1.);
	SECTION("Point angle ctor, horizontal line")
	{
		InfiniteLine l(po, angle);
		REQUIRE(l.isPointCoincident(po, tol));
		REQUIRE(l.isPointCoincident(px, tol));
		REQUIRE(!l.isPointCoincident(py, tol));

		REQUIRE(l.getNormal(po, tol) == py);

		REQUIRE(l.doesLineSegmentIntersect(po, px));
		REQUIRE(l.doesLineSegmentIntersect(po, py));
		REQUIRE(l.doesLineSegmentIntersect(py, px));
		REQUIRE(l.doesLineSegmentIntersect(u, v));
		REQUIRE(!l.doesLineSegmentIntersect(u, py));
	}
	SECTION("y = m*x + b ctor, horizontal line")
	{
		InfiniteLine l(m, b);
		REQUIRE(l.isPointCoincident(po, tol));
		REQUIRE(l.isPointCoincident(px, tol));
		REQUIRE(!l.isPointCoincident(py, tol));

		REQUIRE(l.getNormal(po, tol) == py);

		REQUIRE(l.doesLineSegmentIntersect(po, px));
		REQUIRE(l.doesLineSegmentIntersect(po, py));
		REQUIRE(l.doesLineSegmentIntersect(py, px));
		REQUIRE(l.doesLineSegmentIntersect(u, v));
		REQUIRE(!l.doesLineSegmentIntersect(u, py));
	}
	SECTION("horizontal line ctor")
  {
    InfiniteLine l(0., true);
    REQUIRE(l.isPointCoincident(po, tol));
    REQUIRE(l.isPointCoincident(px, tol));
    REQUIRE(!l.isPointCoincident(py, tol));

    REQUIRE(l.getNormal(po, tol) == py);

    REQUIRE(l.doesLineSegmentIntersect(po, px));
    REQUIRE(l.doesLineSegmentIntersect(po, py));
    REQUIRE(l.doesLineSegmentIntersect(py, px));
    REQUIRE(l.doesLineSegmentIntersect(u, v));
    REQUIRE(!l.doesLineSegmentIntersect(u, py));
  }
	angle = pi/2.;
	SECTION("Point angle ctor, vertical line")
  {
    InfiniteLine l(po, angle);
    REQUIRE(l.isPointCoincident(po, tol));
    REQUIRE(!l.isPointCoincident(px, tol));
    REQUIRE(l.isPointCoincident(py, tol));

    REQUIRE(CGAL::squared_distance(l.getNormal(po, tol), mpx) < tol*tol);

    REQUIRE(l.doesLineSegmentIntersect(po, px));
    REQUIRE(l.doesLineSegmentIntersect(po, py));
    REQUIRE(l.doesLineSegmentIntersect(py, px));
    REQUIRE(l.doesLineSegmentIntersect(u, v));
    REQUIRE(!l.doesLineSegmentIntersect(u, px));
  }
  SECTION("vertical line ctor")
  {
    InfiniteLine l(0., false);
		REQUIRE(l.isPointCoincident(po, tol));
		REQUIRE(!l.isPointCoincident(px, tol));
		REQUIRE(l.isPointCoincident(py, tol));

		REQUIRE(CGAL::squared_distance(l.getNormal(po, tol), mpx) < tol*tol);

		REQUIRE(l.doesLineSegmentIntersect(po, px));
		REQUIRE(l.doesLineSegmentIntersect(po, py));
		REQUIRE(l.doesLineSegmentIntersect(py, px));
		REQUIRE(l.doesLineSegmentIntersect(u, v));
		REQUIRE(!l.doesLineSegmentIntersect(u, px));
  }
}

TEST_CASE("Testing GeometricEntity::LineSegment class", "[GeometricEntity::LineSegment]")
{
  Point_2_base po(0., 0.), px(1.,0.), py(0.,1.), mpx(-1.,0.), mpy(0., -1.);
  double tol = 1e-14;
  Point_2_base u(1., 1.), v(-1., -1.);
	SECTION("2D Line segment")
	{
		LineSegment ls(v, u);
		Point_2_base t(0.5, 0.5);
		REQUIRE(ls.isPointCoincident(po, tol));
		REQUIRE(ls.isPointCoincident(u, tol));
		REQUIRE(ls.isPointCoincident(t, tol));
		REQUIRE(!ls.isPointCoincident(py, tol));

		Point_2_base n(1./sqrt(2), -1./sqrt(2));
		REQUIRE(CGAL::squared_distance(ls.getNormal(po, tol), n) < tol*tol);

		REQUIRE(ls.doesLineSegmentIntersect(po, px));
		REQUIRE(ls.doesLineSegmentIntersect(t, px));
		REQUIRE(ls.doesLineSegmentIntersect(t, py));
		REQUIRE(ls.doesLineSegmentIntersect(t, po));
		REQUIRE(ls.doesLineSegmentIntersect(mpx, mpy));
		REQUIRE(!ls.doesLineSegmentIntersect(mpx, py));
	}
	Point_3_base po3(0.,0.,0.), px3(1.,0.,0.), py3(0.,1.,0.), pz3(0.,0.,1.);
	Point_3_base u3(1.,1.,1.), v3(-1.,-1.,-1.);
	SECTION("3D line segement")
	{
		LineSegment3D ls(v3, u3);
		Point_3_base t(0.5, 0.5, 0.5);
		REQUIRE(ls.isPointCoincident(po3, tol));
		REQUIRE(!ls.isPointCoincident(px3, tol));
		REQUIRE(ls.isPointCoincident(t, tol));
		REQUIRE(ls.isPointCoincident(v3, tol));
		REQUIRE(ls.isPointCoincident(u3, tol));
		REQUIRE(!ls.isPointCoincident(py3, tol));

		REQUIRE(ls.doesLineSegmentIntersect(po3, v3));
		REQUIRE(ls.doesLineSegmentIntersect(po3, u3));
		REQUIRE(ls.doesLineSegmentIntersect(t, px3));
		u3 = Point_3_base(1., 1., 0.);
		v3 = Point_3_base(-1.,-1.,0.);
		REQUIRE(ls.doesLineSegmentIntersect(u3, v3));
		REQUIRE(!ls.doesLineSegmentIntersect(px3, py3));
	}
}

TEST_CASE("Testing GeometricEntity::Polygon2D class","[GeometricEntity::Polygon2D]")
{
	double tol = 1e-14;
	std::vector<Point_2_base> verts = {Point_2_base(-1.,-1.), Point_2_base(1.,-1.), Point_2_base(1.,1.), Point_2_base(-1.,1.)};
	Point_2_base po(0.,0.), t(-1., 0.), u(1.,0.), v(0.,-1.), w(1.,0.);
	Polygon2D p2d(verts);
	SECTION("Coincident points")
	{
		REQUIRE(!p2d.isPointCoincident(po, tol));
		REQUIRE(p2d.isPointCoincident(t, tol));
		REQUIRE(p2d.isPointCoincident(u, tol));
		REQUIRE(p2d.isPointCoincident(v, tol));
		REQUIRE(p2d.isPointCoincident(w, tol));
		for(auto cit = verts.begin(); cit != verts.end(); ++cit)
			REQUIRE(p2d.isPointCoincident(*cit, tol));
	}
	SECTION("Line segements")
	{
		REQUIRE(p2d.doesLineSegmentIntersect(po, t));
		REQUIRE(p2d.doesLineSegmentIntersect(po, u));
		REQUIRE(p2d.doesLineSegmentIntersect(po, v));
		REQUIRE(p2d.doesLineSegmentIntersect(po, w));
	}
	SECTION("Normals")
	{
		REQUIRE(CGAL::squared_distance(p2d.getNormal(t, tol), t) < tol*tol);
		REQUIRE(CGAL::squared_distance(p2d.getNormal(u, tol), u) < tol*tol);
		REQUIRE(CGAL::squared_distance(p2d.getNormal(v, tol), v) < tol*tol);
		REQUIRE(CGAL::squared_distance(p2d.getNormal(w, tol), w) < tol*tol);
	}
}

TEST_CASE("Testing GeometricEntity::InfinitePlane class","[GeometricEntity::InfinitePlane]")
{
	double tol = 1e-14;
	Point_3_base po(0., 0., 0.), n(1., 1., 1.), mn(-1.,-1.,-1.);
	Point_3_base px(1., 0., 0.), py(0., 1., 0.), pz(0., 0., 1.);
	double len = sqrt(3.);
	Point_3_base nhat(1./len, 1./len, 1./len);
	SECTION("Point/normal ctor")
	{
		InfinitePlane ip(po, n);
		// Coincident points
		REQUIRE(ip.isPointCoincident(po, tol));
		REQUIRE(!ip.isPointCoincident(n, tol));
		REQUIRE(ip.isPointCoincident(Point_3_base(-1., 1., 0.), tol));
		// Normal
		REQUIRE(CGAL::squared_distance(ip.getNormal(po, tol), nhat) < tol*tol);
		// Line segment
		REQUIRE(ip.doesLineSegmentIntersect(po, n));
		REQUIRE(ip.doesLineSegmentIntersect(mn, n));
		REQUIRE(!ip.doesLineSegmentIntersect(n, nhat));
	}	
	SECTION("Intercept/plane orientation (XY) ctor")
	{
		InfinitePlane ip(0., GeometryEntityFactory::poXY);
		// Coincident points
		REQUIRE(ip.isPointCoincident(po, tol));
		REQUIRE(!ip.isPointCoincident(n, tol));
		REQUIRE(ip.isPointCoincident(Point_3_base(-1., 1., 0.), tol));
		// Normal
		REQUIRE(CGAL::squared_distance(ip.getNormal(po, tol), pz) < tol*tol);
		// Line segment
		REQUIRE(ip.doesLineSegmentIntersect(po, n));
		REQUIRE(ip.doesLineSegmentIntersect(mn, n));
		REQUIRE(!ip.doesLineSegmentIntersect(n, nhat));
	}
	SECTION("Intercept/plane orientation (YZ) ctor")
  {
		InfinitePlane ip(0., GeometryEntityFactory::poYZ);
		// Coincident points
		REQUIRE(ip.isPointCoincident(po, tol));
		REQUIRE(!ip.isPointCoincident(n, tol));
		REQUIRE(ip.isPointCoincident(Point_3_base(0., 1., 1.), tol));
		// Normal
		REQUIRE(CGAL::squared_distance(ip.getNormal(po, tol), px) < tol*tol);
		// Line segment
		REQUIRE(ip.doesLineSegmentIntersect(po, n));
		REQUIRE(ip.doesLineSegmentIntersect(mn, n));
		REQUIRE(!ip.doesLineSegmentIntersect(n, nhat));
	}
	SECTION("Intercept/plane orientation (ZX) ctor")
	{
		InfinitePlane ip(0., GeometryEntityFactory::poXZ);
		// Coincident points
		REQUIRE(ip.isPointCoincident(po, tol));
		REQUIRE(!ip.isPointCoincident(n, tol));
		REQUIRE(ip.isPointCoincident(Point_3_base(1., 0., 1.), tol));
		// Normal
		REQUIRE(CGAL::squared_distance(ip.getNormal(po, tol), py) < tol*tol);
		// Line segment
		REQUIRE(ip.doesLineSegmentIntersect(po, n));
		REQUIRE(ip.doesLineSegmentIntersect(mn, n));
		REQUIRE(!ip.doesLineSegmentIntersect(n, nhat));	
	}
	SECTION("a*x + b*y + c*z = d ctor")
	{
		InfinitePlane ip(1., 1., 1., 0.);
		// Coincident points
		REQUIRE(ip.isPointCoincident(po, tol));
		REQUIRE(!ip.isPointCoincident(n, tol));
		REQUIRE(ip.isPointCoincident(Point_3_base(-1., 1., 0.), tol));
		// Normal
		REQUIRE(CGAL::squared_distance(ip.getNormal(po, tol), nhat) < tol*tol);
		// Line segment
		REQUIRE(ip.doesLineSegmentIntersect(po, n));
		REQUIRE(ip.doesLineSegmentIntersect(mn, n));
		REQUIRE(!ip.doesLineSegmentIntersect(n, nhat));
	}
	SECTION("Three point ctor")
	{
		InfinitePlane ip(po, Point_3_base(-1., 1., 0.), Point_3_base(0., 0., 1.));
		REQUIRE(ip.isPointCoincident(po, tol));
		REQUIRE(!ip.isPointCoincident(n, tol));
		REQUIRE(ip.isPointCoincident(Point_3_base(0., 0., 2.), tol));
		// Normal
		REQUIRE(CGAL::squared_distance(ip.getNormal(po, tol), Point_3_base(1./sqrt(2), 1./sqrt(2.),0.)) < tol*tol);
		// Line segment
		REQUIRE(ip.doesLineSegmentIntersect(po, n));
		REQUIRE(ip.doesLineSegmentIntersect(mn, n));
		REQUIRE(!ip.doesLineSegmentIntersect(n, nhat));
	}
}

