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

#ifndef GEOMETRICENTITY
#define GEOMETRICENTITY

#include "IO/inputloader.h"
#include "REP/cgal_types.h"
#include <vector>
#include <memory>

//! An enum for defining boundary condition geometry
enum BCType {bctPoint2D, bctVLine, bctHLine, bctLineSeg, bctPolygon2D, bctPoint3D, bctInfinitePlane, bctXYPlane, bctYZPlane, bctXZPlane, bctLineSeg3D, bctUnknown};

//! Light-weight class to handle basic geometric objects such as lines and planes.  Basically a wrapper for CGAL.
/*! This is intended as a helper class for determining which nodes to include for a given boundary condition.
 *  Each derived class contains the object's defining properties, and a function to  determine if a point is
 *  coincident with the object, within some tolerance.
 */
class GeometricEntity
{
public:
	//! Constructor taking the dimension of the geometry and space
	GeometricEntity(unsigned ngd, unsigned nsd) : geomDim(ngd), spatDim(spatDim) {}
	//! Returns whether or not the argument is conincident with this GeometricEntity (2d)
	bool isPointCoincident(const Topologies::Point_2_base& testP) const {return isPointCoincident(testP, 1e-15);}
	//! Returns whether or not the argument is conincident with this GeometricEntity (3d)
	bool isPointCoincident(const Topologies::Point_3_base& testP) const {return isPointCoincident(testP, 1e-15);}
	//! Returns whether or not the argument is conincident with this GeometricEntity within a given tolerance (2d)
	virtual bool isPointCoincident(const Topologies::Point_2_base& testP, const double tol) const = 0;
	//! Returns whether or not the argument is conincident with this GeometricEntity within a given tolerance(3d)
	virtual bool isPointCoincident(const Topologies::Point_3_base& testP, const double tol) const = 0;
	//! Returns a normal vector to this geometry
	//! In 2D, the normal is considered as a cross product with zhat
	virtual Topologies::Point_2_base getNormal(const Topologies::Point_2_base& inP, const double tol) const = 0; 
	//! Returns a normal vector to this geometry
	virtual Topologies::Point_3_base getNormal(const Topologies::Point_3_base& inP, const double tol) const = 0;
	//! Returns whether or not the line segment defined by the arguments intersects wtih this GeometricEntity
	virtual bool doesLineSegmentIntersect(const Topologies::Point_2_base& inP1, const Topologies::Point_2_base& inP2) const = 0;
	//! Returns whether or not the line segment defined by the arguments intersects wtih this GeometricEntity
	virtual bool doesLineSegmentIntersect(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2) const = 0;
	//! Returns a new copy of this GeometricEntity
	virtual std::unique_ptr<GeometricEntity> clone() const = 0;
	virtual ~GeometricEntity() {}
	
	unsigned getGeomDim() const {return geomDim;}
	unsigned getSpatDim() const {return spatDim;}
private:
	unsigned geomDim, spatDim;
};

//! Namespace containing enums and functions to generate GeometricEntity objects
namespace GeometryEntityFactory
{
	//! Enum defining 3d infinite plane orientation
	enum PlaneOrientation {poXY, poYZ, poXZ};
	//! Factory for creating GeometricEntity objects
	std::unique_ptr<GeometricEntity> createGeometricEntity(const BCType inBCT, const pugi::xml_node& rootNode);
	//! Factory for creating GeometricEntity objects
	std::unique_ptr<GeometricEntity> createGeometricEntity(const BCType inBCT, const pugi::xml_node& rootNode, 
																													const Topologies::DiscretizationParameters& discParams);
	//! Returns a 2d point from a 3d point using coordinates as defined by the plane orientation
	Topologies::Point_2_base point2Dfrom3D(const Topologies::Point_3_base& inPt3, GeometryEntityFactory::PlaneOrientation projPO);
}

//! Point in 3D space (used for both 2d and 3d)
class Point : public GeometricEntity
{
public:
	//! Constructor taking a Point_2_base
	Point(const Topologies::Point_2_base& inP);
	//! Constructor taking a Point_3_base
	Point(const Topologies::Point_3_base& inP);
	virtual bool isPointCoincident(const Topologies::Point_2_base& testP, const double tol) const;
	virtual bool isPointCoincident(const Topologies::Point_3_base& testP, const double tol) const;
	virtual Topologies::Point_2_base getNormal(const Topologies::Point_2_base& inP, const double tol) const {return Topologies::Point_2_base(0.,0.);}
	virtual Topologies::Point_3_base getNormal(const Topologies::Point_3_base& inP, const double tol) const {return Topologies::Point_3_base(0.,0.,0.);}
	virtual bool doesLineSegmentIntersect(const Topologies::Point_2_base& inP1, const Topologies::Point_2_base& inP2) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2) const;
	virtual std::unique_ptr<GeometricEntity> clone() const {return std::unique_ptr<GeometricEntity>(new Point(*this));}
private:
	Topologies::Point_3_base p;
};

//! Infinite line in 2d space: Defined by a point and an angle
class InfiniteLine : public GeometricEntity
{
public:
	//! Constructor defining an infinite line from a point and an angle
	InfiniteLine(const Topologies::Point_2_base& inOrigin, const double inAngle);
	//! Constructor defining an infinite line from a slob and an offset
	InfiniteLine(const double m, const double b);
	//! Constructor defining a horizontal or verticle infinite line 
	InfiniteLine(const double intercept, bool isHorizontal);
	virtual bool isPointCoincident(const Topologies::Point_2_base& testP, const double tol) const;
	virtual Topologies::Point_2_base getNormal(const Topologies::Point_2_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_2_base& inP1, const Topologies::Point_2_base& inP2) const;
	virtual std::unique_ptr<GeometricEntity> clone() const {return std::unique_ptr<GeometricEntity>(new InfiniteLine(*this));}
private:
	Topologies::Point_2_base origin;
	double angle;
private:
	// Hidden, no implementation
	virtual bool isPointCoincident(const Topologies::Point_3_base& testP, const double tol) const;
	virtual Topologies::Point_3_base getNormal(const Topologies::Point_3_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2) const;
};

//! Line segment in 2d space: Defined by two end points
class LineSegment : public GeometricEntity
{
public:
	//! Constructor defining a line segment from two points
	LineSegment(const Topologies::Point_2_base& inP1, const Topologies::Point_2_base& inP2);
	virtual bool isPointCoincident(const Topologies::Point_2_base& testP, const double tol) const;
	virtual Topologies::Point_2_base getNormal(const Topologies::Point_2_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_2_base& inP1, const Topologies::Point_2_base& inP2) const;
	virtual std::unique_ptr<GeometricEntity> clone() const {return std::unique_ptr<GeometricEntity>(new LineSegment(*this));}

	std::vector<Topologies::Point_2_base> getDefiningPoints() const;
private:
	Topologies::Point_2_base p1, p2;
private:
	// Hidden, no implementation
	virtual bool isPointCoincident(const Topologies::Point_3_base& testP, const double tol) const;
	virtual Topologies::Point_3_base getNormal(const Topologies::Point_3_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2) const;
};

//! Polygon in 2d space: Polygon defined by a set of points, applies only to the boundary
class Polygon2D : public GeometricEntity
{
public:
	//! Constructor defining a polygon from a vector of vertices (must be in order!)
	Polygon2D(const std::vector<Topologies::Point_2_base>& inVerts);
	~Polygon2D(){}
	Polygon2D(const Polygon2D& rhsP);
	Polygon2D(Polygon2D&& other);
	Polygon2D& operator=(Polygon2D rhsP);

	virtual bool isPointCoincident(const Topologies::Point_2_base& testP, const double tol) const;
	virtual Topologies::Point_2_base getNormal(const Topologies::Point_2_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_2_base& inP1, const Topologies::Point_2_base& inP2) const;
	virtual std::unique_ptr<GeometricEntity> clone() const {return std::unique_ptr<GeometricEntity>(new Polygon2D(*this));}

	std::vector<Topologies::Point_2_base> getDefiningPoints() const;
private:
	friend void swap(Polygon2D& first, Polygon2D& second) { std::swap(first.lineSegVec, second.lineSegVec); }

	std::vector<std::unique_ptr<LineSegment> > lineSegVec;
private:
	// Hidden, no implementation
	virtual bool isPointCoincident(const Topologies::Point_3_base& testP, const double tol) const;
	virtual Topologies::Point_3_base getNormal(const Topologies::Point_3_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2) const;
};

//! Line segment in 2d space: Defined by two end points
class LineSegment3D : public GeometricEntity
{
public:
	//! Constructor defining a line segment from two points
	LineSegment3D(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2);
	virtual bool isPointCoincident(const Topologies::Point_3_base& testP, const double tol) const;
	virtual Topologies::Point_3_base getNormal(const Topologies::Point_3_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2) const;
	virtual std::unique_ptr<GeometricEntity> clone() const {return std::unique_ptr<GeometricEntity>(new LineSegment3D(*this));}

private:
	Topologies::Point_3_base p1, p2;
private:
	// Hidden, no implementation
	virtual bool isPointCoincident(const Topologies::Point_2_base& testP, const double tol) const;
	virtual Topologies::Point_2_base getNormal(const Topologies::Point_2_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_2_base& inP1, const Topologies::Point_2_base& inP2) const;
};

//! Infinite plane in 3d space: Defined by a point and a normal
class InfinitePlane : public GeometricEntity
{
public:
	//! Constructor defining an infinite plane from a point and a normal
	InfinitePlane(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inN1);
	//! Constructor defining an infinite plane in the x-y, y-z, or z-x planes
	InfinitePlane(const double intercept, const GeometryEntityFactory::PlaneOrientation inPO, double normalDir = 1.);
	//! Constructor defining an infinite plane from its defining equation a*x + b*y + c*z = d
	InfinitePlane(const double a, const double b, const double c, const double d);
	//! Constructor defining an infinite plane from three (non-collinear) points
	InfinitePlane(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2, const Topologies::Point_3_base& inP3);
	virtual bool isPointCoincident(const Topologies::Point_3_base& testP, const double tol) const;
	virtual Topologies::Point_3_base getNormal(const Topologies::Point_3_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_3_base& inP1, const Topologies::Point_3_base& inP2) const;
	virtual std::unique_ptr<GeometricEntity> clone() const {return std::unique_ptr<GeometricEntity>(new InfinitePlane(*this));}

	Topologies::Point_3_base getIntercept() const {return p1;};
private:
	Topologies::Point_3_base p1;
	Topologies::Mesh_K::Vector_3 n1;
private:
	// Hidden
	virtual bool isPointCoincident(const Topologies::Point_2_base& testP, const double tol) const;
	virtual Topologies::Point_2_base getNormal(const Topologies::Point_2_base& inP, const double tol) const;
	virtual bool doesLineSegmentIntersect(const Topologies::Point_2_base& inP1, const Topologies::Point_2_base& inP2) const;
};

#endif
