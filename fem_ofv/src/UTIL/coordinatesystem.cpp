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

#include "coordinatesystem.h"
#include "point2d.h"
#include "point3d.h"

namespace CoordinateSystem{
namespace{
	template<typename PointType>
	double phi(PointType const& p)
	{
		return atan2(p.y, p.x);
	}
	double theta(Point3D const& p)
	{
		double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
		return acos(p.z/r);
	}
}

  template<typename PointType>
  PointType cylToCart(PointType const& p)
	{
		double r = p.x;
		double ph = p.y;
		PointType q = p;  // Copies z for 3d
		q.x = r*cos(ph);
		q.y = r*sin(ph);
		return q;
	}

  template<typename PointType>
  PointType cylVectorToCartVector(PointType const& v, PointType const& p)
	{
		double ph = phi(p);
		PointType q = v; // Copies z coordinate for 3d
		q.x = v.x*cos(ph) - v.y*sin(ph);
		q.y = v.x*sin(ph) + v.y*cos(ph);
		return q;
	}

  Point3D sphrToCart(Point3D const& p)
	{
		double r = p.x;
		double th = p.y;
		double ph = p.z;
		return Point3D(r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th));
	}

  Point3D sphrVectorToCartVector(Point3D const& v, Point3D const& p)
	{
		double ph = phi(p);
		double th = theta(p);
		return Point3D(/*x=*/ v.x*sin(th)*cos(ph) + v.y*cos(th)*cos(ph) - v.z*sin(ph),
									 /*y=*/ v.x*sin(th)*sin(ph) + v.y*cos(th)*sin(ph) + v.z*cos(ph),
									 /*z=*/ v.x*cos(th) - v.y*sin(th));
	}

  Point3D convertVector(Point3D const& v, Point3D const& p, Type ct)
	{
		if(ct == Type::spherical)
			return sphrVectorToCartVector(v, p);
		else if(ct == Type::cylindrical)
			return cylVectorToCartVector(v, p);
		return v;
	}
}

template Point2D CoordinateSystem::cylVectorToCartVector<Point2D>(Point2D const& v, Point2D const& p);
template Point3D CoordinateSystem::cylVectorToCartVector<Point3D>(Point3D const& v, Point3D const& p);
template Point2D CoordinateSystem::cylToCart<Point2D>(Point2D const& p);
template Point3D CoordinateSystem::cylToCart<Point3D>(Point3D const& p);


