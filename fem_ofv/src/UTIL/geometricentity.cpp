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

#include "geometricentity.h"
#include <memory>

const double defaultTol = 1e-15;

using namespace Topologies;

// ******** Point in 3D space (used for both 2d and 3d)

Point::Point(const Point_2_base& inP) : 
	GeometricEntity(0, 2),
	p(Point_3_base(inP.x(), inP.y(), 0.))
{
}

Point::Point(const Point_3_base& inP) :
	GeometricEntity(0, 3),
	p(inP) 
{
}

bool Point::isPointCoincident(const Point_2_base& testP, const double tol) const
{
	return isPointCoincident(Point_3_base(testP.x(), testP.y(), 0.), tol);
}

bool Point::isPointCoincident(const Point_3_base& testP, const double tol) const
{
	return CGAL::squared_distance(testP, p) < tol*tol;
}

bool Point::doesLineSegmentIntersect(const Point_2_base& inP1, const Point_2_base& inP2) const
{
	Mesh_K::Segment_2 tmpSeg2(inP1, inP2);
	Point_2_base tmpp(p.x(), p.y());
	return CGAL::do_intersect(tmpSeg2, tmpp);
}

bool Point::doesLineSegmentIntersect(const Point_3_base& inP1, const Point_3_base& inP2) const
{
	const Mesh_K::Segment_3 tmpSeg3(inP1, inP2);
	return tmpSeg3.has_on(p);
}

// ******** Infinite line in 2d space: Defined by a point and an angle

InfiniteLine::InfiniteLine(const Point_2_base& inOrigin, const double inAngle) : 
	GeometricEntity(1, 2),
	origin(inOrigin),
	angle(inAngle)
{
} // Point/angle form

InfiniteLine::InfiniteLine(const double m, const double b) : // Standard form, y = m*x + b
	GeometricEntity(1, 2),
	origin(Point_2_base(0., b)),
	angle(atan(m))
{
}

InfiniteLine::InfiniteLine(const double intercept, bool isHorizontal) :  // Horizontal/vertical line form
	GeometricEntity(1, 2),
	origin(isHorizontal ? Point_2_base(0., intercept) : Point_2_base(intercept, 0.)), 
	angle(isHorizontal ? 0. : PI*0.5)
{
}

bool InfiniteLine::isPointCoincident(const Point_2_base& testP, const double tol) const
{
	Mesh_K::Vector_2 tmpVec(cos(angle), sin(angle));
	Mesh_K::Line_2 tmpLine(origin, tmpVec);
	return CGAL::squared_distance(testP, tmpLine) < tol*tol;
}

bool InfiniteLine::isPointCoincident(const Point_3_base& testP, const double tol) const
{
	return isPointCoincident(Point_2_base(testP.x(), testP.y()), tol);
}

Point_2_base InfiniteLine::getNormal(const Point_2_base& inP, const double tol) const
{
	return Point_2_base(-sin(angle), cos(angle));
}

Point_3_base InfiniteLine::getNormal(const Point_3_base& inP, const double tol) const
{
	Point_2_base tmp = getNormal(Point_2_base(inP.x(), inP.y()), tol);
	return Point_3_base(tmp.x(), tmp.y(), 0.);
}

bool InfiniteLine::doesLineSegmentIntersect(const Point_2_base& inP1, const Point_2_base& inP2) const
{
	Mesh_K::Segment_2 tmpSeg2(inP1, inP2);
	Mesh_K::Vector_2 tmpVec(cos(angle), sin(angle));
  Mesh_K::Line_2 tmpLine(origin, tmpVec);
  return CGAL::do_intersect(tmpSeg2, tmpLine);
}

bool InfiniteLine::doesLineSegmentIntersect(const Point_3_base& inP1, const Point_3_base& inP2) const
{
	// Not implemented
	return false;
}

// ********** Line segment in 2d space: Defined by two end points

LineSegment::LineSegment(const Point_2_base& inP1, const Point_2_base& inP2) :
	GeometricEntity(1, 2),
	p1(inP1),
	p2(inP2)
{
}

bool LineSegment::isPointCoincident(const Point_2_base& testP, const double tol) const
{
	Mesh_K::Segment_2 tmpSeg2(p1, p2);
	return CGAL::squared_distance(tmpSeg2, testP) < tol*tol;
}

bool LineSegment::isPointCoincident(const Point_3_base& testP, const double tol) const
{
	return isPointCoincident(Point_2_base(testP.x(), testP.y()), tol);
}

Point_2_base LineSegment::getNormal(const Point_2_base& inP, const double tol) const
{
	Mesh_K::Vector_2 v1 = p2 - p1;
	Mesh_K::Vector_2 n1(v1.y(), -v1.x());
	return CGAL::ORIGIN + n1/sqrt(n1.squared_length());
}

Point_3_base LineSegment::getNormal(const Point_3_base& inP, const double tol) const
{
	Point_2_base tmp = getNormal(Point_2_base(inP.x(), inP.y()), tol);
	return Point_3_base(tmp.x(), tmp.y(), 0.);
}

bool LineSegment::doesLineSegmentIntersect(const Point_2_base& inP1, const Point_2_base& inP2) const
{
	Mesh_K::Segment_2 tmpSeg1(p1, p2);
	Mesh_K::Segment_2 tmpSeg2(inP1, inP2);
	return CGAL::do_intersect(tmpSeg1, tmpSeg2);
}

bool LineSegment::doesLineSegmentIntersect(const Point_3_base& inP1, const Point_3_base& inP2) const
{
	// Not implemented
	return false;
}

std::vector<Point_2_base> LineSegment::getDefiningPoints() const
{
	std::vector<Point_2_base> ptVec(2);
	ptVec[0] = p1;
	ptVec[1] = p2;
	return ptVec;
}

// ********** Polygon in 2d space: Polygon defined by a set of points, applies only to the boundary

Polygon2D::Polygon2D(const std::vector<Point_2_base>& inVerts) :
	GeometricEntity(2, 2)
{
	assert(inVerts.size() >= 3);
	for(unsigned k = 0; k < inVerts.size(); ++k)
	{
		unsigned kp1 = (k + 1) % inVerts.size();
		lineSegVec.push_back(std::unique_ptr<LineSegment>(new LineSegment(inVerts[k], inVerts[kp1])));
	}
}

Polygon2D::Polygon2D(const Polygon2D& rhsP) :
	GeometricEntity(rhsP)
{
	for(unsigned k = 0; k < rhsP.lineSegVec.size(); ++k)
		lineSegVec.push_back(std::unique_ptr<LineSegment>(new LineSegment(*rhsP.lineSegVec[k])));
}

Polygon2D::Polygon2D(Polygon2D&& other) :
	GeometricEntity(other)
{
	swap(*this, other);
}

Polygon2D& Polygon2D::operator=(Polygon2D rhsP)
{
	swap(*this, rhsP);
	return *this;
}

bool Polygon2D::isPointCoincident(const Point_2_base& testP, const double tol) const
{
	// Use class LineSegment to do checks
	if(lineSegVec.size() < 3)
		return false;
	for(unsigned k = 0; k < lineSegVec.size(); ++k)
	{
		if(lineSegVec[k]->isPointCoincident(testP, tol))
			return true;
	}
	return false;
}

bool Polygon2D::isPointCoincident(const Point_3_base& testP, const double tol) const
{
	return isPointCoincident(Point_2_base(testP.x(), testP.y()), tol);
}

Point_2_base Polygon2D::getNormal(const Point_2_base& inP, const double tol) const
{
	// Loops over all sides and takes sum of normal if point is coincident w/ 2 sides (i.e. a corner)
	Mesh_K::Vector_2 nhat(0.,0.);
	for(unsigned k = 0; k < lineSegVec.size(); ++k)
	{
		if(lineSegVec[k]->isPointCoincident(inP, tol))
			nhat = nhat + (lineSegVec[k]->getNormal(inP, tol) - CGAL::ORIGIN);
	}
	nhat = nhat/sqrt(nhat.squared_length());
	return CGAL::ORIGIN + nhat;
}

Point_3_base Polygon2D::getNormal(const Point_3_base& inP, const double tol) const
{
	Point_2_base tmp = getNormal(Point_2_base(inP.x(), inP.y()), tol);
	return Point_3_base(tmp.x(), tmp.y(), 0.);
}

bool Polygon2D::doesLineSegmentIntersect(const Point_2_base& inP1, const Point_2_base& inP2) const
{
	// Use class LineSegment to do checks
	if(lineSegVec.size() < 3)
		return false;
	for(unsigned k = 0; k < lineSegVec.size(); ++k)
	{
		if(lineSegVec[k]->doesLineSegmentIntersect(inP1, inP2))
			return true;
	}
	return false;
}

bool Polygon2D::doesLineSegmentIntersect(const Point_3_base& inP1, const Point_3_base& inP2) const
{
	// Not implemented
	return false;
}

std::vector<Point_2_base> Polygon2D::getDefiningPoints() const
{
	std::vector<Point_2_base> ptVec(lineSegVec.size());
	for(unsigned k = 0; k < lineSegVec.size(); ++k)
	{
		std::vector<Point_2_base> tmppts = lineSegVec[k]->getDefiningPoints();
		ptVec[k] = tmppts[0];
	}
	return ptVec;
}

// ********* Line segment in 3d space

LineSegment3D::LineSegment3D(const Point_3_base& inP1, const Point_3_base& inP2):
	GeometricEntity(1, 3),
	p1(inP1),
	p2(inP2)
{
}

bool LineSegment3D::isPointCoincident(const Point_3_base& testP, const double tol) const
{
	Mesh_K::Segment_3 tmpSeg2(p1, p2);
  return CGAL::squared_distance(tmpSeg2, testP) < tol*tol;
}

Point_3_base LineSegment3D::getNormal(const Point_3_base& inP, const double tol) const
{
	// This is not correct, though nothing in this code uses getNormal
	Point_2_base p12(p1.x(), p1.y()), p22(p2.x(), p2.y());
	Mesh_K::Vector_2 v1 = p22 - p12;
	Mesh_K::Vector_2 n1(v1.y(), -v1.x());
	Point_2_base tmp = CGAL::ORIGIN + n1/sqrt(n1.squared_length());
	return Point_3_base(tmp.x(), tmp.y(), 0.);
}

bool LineSegment3D::doesLineSegmentIntersect(const Point_3_base& inP1, const Point_3_base& inP2) const
{
	Mesh_K::Segment_3 tmpSeg1(p1, p2);
	Mesh_K::Segment_3 tmpSeg2(inP1, inP2);
	return CGAL::do_intersect(tmpSeg1, tmpSeg2);
}

bool LineSegment3D::isPointCoincident(const Point_2_base& testP, const double tol) const
{
	return false;
}

Point_2_base LineSegment3D::getNormal(const Point_2_base& inP, const double tol) const
{
	return Point_2_base(0., 0.);
}

bool LineSegment3D::doesLineSegmentIntersect(const Point_2_base& inP1, const Point_2_base& inP2) const
{
	return false;
}

// ********* Infinite plane in 3d space: Defined by a point and a normal

InfinitePlane::InfinitePlane(const Point_3_base& inP1, const Point_3_base& inN1) : 
	GeometricEntity(2, 3),
	p1(inP1), 
	n1(CGAL::ORIGIN, inN1)
{
}

InfinitePlane::InfinitePlane(const double intercept, const GeometryEntityFactory::PlaneOrientation inPO, double normalDir):
	GeometricEntity(2, 3),
	p1(inPO == GeometryEntityFactory::poXY ? Point_3_base(0., 0., intercept) : (inPO == GeometryEntityFactory::poYZ ? Point_3_base(intercept, 0., 0.) : Point_3_base(0., intercept, 0.))),
	n1(inPO == GeometryEntityFactory::poXY ? Mesh_K::Vector_3(0., 0., 1.) : (inPO == GeometryEntityFactory::poYZ ? Mesh_K::Vector_3(1., 0., 0.) : Mesh_K::Vector_3(0., 1., 0.)))
{
	n1 = (normalDir >= 0. ? 1. : -1.)*n1;
}

InfinitePlane::InfinitePlane(const double a, const double b, const double c, const double d):
	GeometricEntity(2, 3),
	p1(c != 0. ? Point_3_base(0., 0., d/c) : (b != 0. ? Point_3_base(0., d/b, 0.) : Point_3_base(d/a, 0., 0.))),
	n1(a, b, c)
{
}

InfinitePlane::InfinitePlane(const Point_3_base& inP1, const Point_3_base& inP2, const Point_3_base& inP3):
	GeometricEntity(2, 3),
	p1(inP1),
	n1(CGAL::cross_product(inP2 - inP1, inP3 - inP1))
{
}

bool InfinitePlane::isPointCoincident(const Point_2_base& testP, const double tol) const 
{ 
	return isPointCoincident(Point_3_base(testP.x(), testP.y(), 0.), tol); 
}

bool InfinitePlane::isPointCoincident(const Point_3_base& testP, const double tol) const
{
	Mesh_K::Plane_3 tmpPlane(p1,n1);
	return CGAL::squared_distance(testP, tmpPlane) < tol*tol;
}

Point_2_base InfinitePlane::getNormal(const Point_2_base& inP, const double tol) const
{
	Point_3_base tmp = getNormal(Point_3_base(inP.x(), inP.y(), 0.), tol);
	return Point_2_base(tmp.x(),tmp.y());
}

Point_3_base InfinitePlane::getNormal(const Point_3_base& inP, const double tol) const
{
	double len = sqrt(n1.squared_length());
	return Point_3_base(n1.x()/len, n1.y()/len, n1.z()/len);
}

bool InfinitePlane::doesLineSegmentIntersect(const Point_2_base& inP1, const Point_2_base& inP2) const
{
	// Not implemented
	return false;
}

bool InfinitePlane::doesLineSegmentIntersect(const Point_3_base& inP1, const Point_3_base& inP2) const
{
	Mesh_K::Plane_3 tmpPlane(p1,n1);
	Mesh_K::Segment_3 tmpSeg3(inP1, inP2);
	return CGAL::do_intersect(tmpPlane, tmpSeg3);
}

// Factory function to generate GeometricEntity objects from BCType and appropriate text inputs from Parser
/* All constructors:
	Point(const Point_2_base& inP)
	Point(const Point_3_base& inP)
	InfiniteLine(const Point_2_base& inOrigin, const double inAngle)
	InfiniteLine(const double m, const double b)
	InfiniteLine(const double intercept, bool isHorizontal)
	LineSegment(const Point_2_base& inP1, const Point_2_base& inP2)
	InfinitePlane(const Point_3_base& inP1, const Point_3_base& inN1)
	InfinitePlane(const double intercept, const PlaneOrientation inPO)
	InfinitePlane(const double a, const double b, const double c, const double d)
	Polygon(const std::vector<Point_3_base>& inVerts)
*/
namespace GeometryEntityFactory
{
	namespace
	{
		using namespace InputLoader;

		Point_2_base readPoint2(const pugi::xml_node& rootNode)
		{
			double x = readDoubleAttribute(rootNode, "x");
      double y = readDoubleAttribute(rootNode, "y");
			return Point_2_base(x, y);
		}
		Point_3_base readPoint3(const pugi::xml_node& rootNode)
    {
			double x = readDoubleAttribute(rootNode, "x");
			double y = readDoubleAttribute(rootNode, "y");
			double z = readDoubleAttribute(rootNode, "z");
			return Point_3_base(x, y, z);
    }
		std::pair<Point_2_base, Point_2_base> readLineSegment2(const pugi::xml_node& rootNode)
		{
			pugi::xml_node cn = rootNode.child("point_2");
			std::pair<Point_2_base, Point_2_base> outLS;
			outLS.first = readPoint2(cn);
			cn = cn.next_sibling("point_2");
			outLS.second = readPoint2(cn);
			return outLS;
		}
		std::pair<Point_3_base, Point_3_base> readLineSegment3(const pugi::xml_node& rootNode)
    {
			pugi::xml_node cn = rootNode.child("point_3");
			std::pair<Point_3_base, Point_3_base> outLS;
			outLS.first = readPoint3(cn);
			cn = cn.next_sibling("point_3");
			outLS.second = readPoint3(cn);
			return outLS;
    }
		std::vector<Point_2_base> readPolygon2(const pugi::xml_node& rootNode)
		{
			std::vector<Point_2_base> ptVec;
			for(pugi::xml_node cn = rootNode.child("point_2"); cn; cn = cn.next_sibling("point_2"))
				ptVec.push_back(readPoint2(cn));
			return ptVec;
  	}
		std::vector<Point_3_base> readPolygon3(const pugi::xml_node& rootNode)
		{
			std::vector<Point_3_base> ptVec;
			for(pugi::xml_node cn = rootNode.child("point_3"); cn; cn = cn.next_sibling("point_3"))
				ptVec.push_back(readPoint3(cn));
			return ptVec;
		}
	}

std::unique_ptr<GeometricEntity> createGeometricEntity(const BCType inBCT, const pugi::xml_node& rootNode, 
																												const DiscretizationParameters& discParams)
{
	using namespace InputLoader;

	std::unique_ptr<GeometricEntity> outGE;
	if(inBCT == bctPoint2D)
	{
		Point_2_base p = readPoint2(rootNode.child("point_2"));
		outGE = std::unique_ptr<GeometricEntity>(new Point(p));
	}
	else if(inBCT == bctVLine || inBCT == bctHLine)
	{
		double intercept = readDoublePCData(rootNode.child("intercept"));
		outGE = std::unique_ptr<GeometricEntity>(new InfiniteLine(intercept, inBCT == bctHLine));
	}
	else if(inBCT == bctLineSeg)
	{
		std::pair<Point_2_base, Point_2_base> ls = readLineSegment2(rootNode);
		outGE = std::unique_ptr<GeometricEntity>(new LineSegment(ls.first, ls.second));
	}
	else if(inBCT == bctLineSeg3D)
	{
		std::pair<Point_3_base, Point_3_base> ls = readLineSegment3(rootNode);
		outGE = std::unique_ptr<GeometricEntity>(new LineSegment3D(ls.first, ls.second));
	}
	else if(inBCT == bctPolygon2D)
	{
		std::vector<Point_2_base> vertices = readPolygon2(rootNode);
		outGE = std::unique_ptr<GeometricEntity>(new Polygon2D(vertices));
	}
	else if(inBCT == bctPoint3D)
	{
		Point_3_base p = readPoint3(rootNode.child("point_3"));
		outGE = std::unique_ptr<GeometricEntity>(new Point(p));
	}
	else if(inBCT == bctXYPlane || inBCT == bctYZPlane || inBCT == bctXZPlane)
	{
		double intercept = readDoublePCData(rootNode.child("intercept"));
		PlaneOrientation curPO = inBCT == bctXYPlane ? poXY : (inBCT == bctYZPlane ? poYZ : poXZ);
		outGE = std::unique_ptr<GeometricEntity>(new InfinitePlane(intercept, curPO));
	}
	else if(inBCT == bctInfinitePlane)
	{
		Point_3_base p = readPoint3(rootNode.child("point"));
		Point_3_base n = readPoint3(rootNode.child("normal"));;
		outGE = std::unique_ptr<GeometricEntity>(new InfinitePlane(p, n));
	}
	return outGE;
}

std::unique_ptr<GeometricEntity> createGeometricEntity(const BCType inBCT, const pugi::xml_node& rootNode)
{
	return createGeometricEntity(inBCT, rootNode, DiscretizationParameters());
}

Point_2_base point2Dfrom3D(const Point_3_base& inPt3, GeometryEntityFactory::PlaneOrientation projPO)
{
	double val1 = 0., val2 = 0.;
	if(projPO == GeometryEntityFactory::poXY)
	{
		val1 = inPt3.x();
		val2 = inPt3.y();
	}
	else if(projPO == GeometryEntityFactory::poYZ)
	{
		val1 = inPt3.y();
		val2 = inPt3.z();
	}
	else if(projPO == GeometryEntityFactory::poXZ)
	{
		val1 = inPt3.x();
		val2 = inPt3.z();
	}
	return Point_2_base(val1, val2);
}
}
