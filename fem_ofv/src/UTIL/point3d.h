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

#ifndef POINT3D_H
#define POINT3D_H

#include "UTIL/topologiesdefs.h"
#include "point2d.h"
#include <iostream>
#include <assert.h>

//! A basic implementation of a point in 3d space
/*! The class has publically accesible coordinates (x,y,z) and numerous functions
 *  to perform operations on 3d points.  
 */
struct Point3D
{
public:
	//! Point coordinates
	double x, y, z;
	//! Default constructor initializes all coordinates to 0
	Point3D();
	//! Constructor sets the x and y coordinates to the given arguments, and z to 0
	Point3D(double inX, double inY);
	//! Constructor that sets all coordinates to the given arguments
	Point3D(double inX, double inY, double inZ);
	Point3D(const Point3D& point); 
	Point3D(const Point3D* point);
	//! Converts a CGAL Point_3 to a Point3D
	Point3D(const Topologies::Point_3_base& point);
	Point3D& operator=(const Point3D& rhs);
	//! Returns the negation of this point, i.e. (-x, -y, -z)
	Point3D operator-() const;
	//! Returns the sum of this point and the argument p
	Point3D operator+(const Point3D& p) const;
	//! Returns the difference of this point and the argument p
	Point3D operator-(const Point3D& p) const;
	//! Sums the argument and stores result in this Point3D
	const Point3D& operator+=(const Point3D& p);
	//! Computes the difference of the argument and stores result in this Point3D
	const Point3D& operator-=(const Point3D& p);
	//! Scalar multiplication of a Point3D
	friend Point3D operator*(double a, const Point3D& p);
	//! Scalar multiplication of a Point3D
	friend Point3D operator*(const Point3D& p, double a); 
	//! Scalar product (dot product) of two Point3D objects
	double operator*(const Point3D& p) const;
	//! Returns the cross product with the argument p
	Point3D crossProduct(const Point3D& p) const;
	//! Returns the cross product of p1 and p2
	friend Point3D crossProduct(const Point3D& p1, const Point3D& p2);
	//! Returns the area of the triangle defined by the three arguments
	friend double triangleArea(const Point3D& p1, const Point3D& p2, const Point3D& p3);
	//! Returns the volume of the tetrahedron formed by the four arguments
	friend double tetrahedronVolume(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4);
	//! Returns the centroid of the tetrahedron formed by the four arguments
	friend Point3D centroid(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4);
	//! Returns the solid angle of the 3 vectors formed by a-o, b-o, c-o
	friend double solidAngle(const Point3D& o, const Point3D& a, const Point3D& b, const Point3D& c);
	//! Computes and stores the result of scalar multiplication with argument a
	const Point3D& operator*=(double a);
	//! Scalar division
	Point3D operator/(double a) const;
	//! Computes and stores the result of scalar division with argument a
	const Point3D& operator/=(double a);
	//! Returns the length of the vector connecting argument v and the origin
	friend double abs(const Point3D& v);
	//! Checks lexographical order of this Point3D and the argument.
	/*! Ordering is given by first checking x-coordinate.  If the x-coordinates are equal, the y-coordinate is used, and so on.*/
	bool operator<(const Point3D& point) const;
	//! Checks equality of the point's coordinates
	bool operator==(const Point3D& point) const;
	//! Checks lexographical order of this Point3D and the argument.
	bool operator<=(const Point3D& point) const;
	//! Checks lexographical order of this Point3D and the argument.
	bool operator>=(const Point3D& point) const;
	//! Checks inequality of the point's coordinates
	bool operator!=(const Point3D& point) const;
	//! Checks lexographical order of this Point3D and the argument.
	bool operator>(const Point3D& point) const;
	//! Text output to the specified output stream.  Outputs the coordinate values
	friend std::ostream& operator<<(std::ostream& theStream, Point3D thePoint);
	//! Text input from the specified input stream.  Sets all three coordinates
	friend std::istream& operator>>(std::istream& theStream, Point3D& thePoint);
	//! Computes and returns the vector norm of the cross product of this Point3D and the argument
	double crossProductNorm(const Point3D& point) const;
	//! Returns the dimension of the space this object lives in
	unsigned dim() const {return 3;}
	//! Coordinate access (const)
	double operator[](unsigned k) const;
private:
	static double tolx, toly, tolz;
};

inline
double Point3D::operator[](unsigned k) const
{
	if(k == 0) 
		return x;
	else if(k == 1)
		return y;
	else if(k == 2)
		return z;
	return 0.;
}

inline
Point3D::Point3D():
	x(0.), y(0.), z(0.)
{
}

inline
Point3D::Point3D(double inX, double inY):
	x(inX), y(inY), z(0.)
{
}

inline
Point3D::Point3D(double inX, double inY, double inZ):
	x(inX), y(inY), z(inZ)
{
}

inline 
Point3D::Point3D(const Point3D& point):
	x(point.x), y(point.y), z(point.z)
{
}
 
inline 
Point3D::Point3D(const Point3D* point):
x(point->x), y(point->y), z(point->z)
{
}

inline
Point3D::Point3D(const Topologies::Point_3_base& point) :
  x(point.x()), y(point.y()), z(point.z())
{
}

inline
Point3D Point3D::operator-() const
{
	return Point3D(-x, -y, -z);
}

inline
Point3D Point3D::operator+(const Point3D& p) const
{
	return Point3D(x + p.x, y + p.y, z + p.z);
}

inline
Point3D Point3D::operator-(const Point3D& p) const
{
	return Point3D(x - p.x, y - p.y, z - p.z);
}

inline 
const Point3D& Point3D::operator+=(const Point3D& p)
{
	x += p.x;
	y += p.y;
	z += p.z;
	return *this;
}

inline 
const Point3D& Point3D::operator-=(const Point3D& p)
{
	x -= p.x;
	y -= p.y;
	z -= p.z;
	return *this;
}

inline
Point3D operator*(double a, const Point3D& p)
{
	return Point3D(a*p.x, a*p.y, a*p.z);
}

inline
Point3D operator*(const Point3D& p, double a)
{
	return Point3D(a*p.x, a*p.y, a*p.z);
}

inline 
const Point3D& Point3D::operator*=(double a)
{
	x *= a;
	y *= a;
	z *= a;
	return *this;
}

inline
Point3D Point3D::operator/(double a) const
{
	return Point3D(x/a, y/a, z/a);
}

inline 
const Point3D& Point3D::operator/=(double a)
{
	x /= a;
	y /= a;
	z /= a;
	return *this;
}

inline 
double abs(const Point3D& v)
{
	return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline 
bool Point3D::operator==(const Point3D& point) const
{
	return (TOL_EQ(x, point.x, tolx)) 
		  && (TOL_EQ(y, point.y, toly))
			&& (TOL_EQ(z, point.z, tolz));
}

inline 
bool Point3D::operator<=(const Point3D& point) const
{
	return (*this < point) || (*this == point);
}

inline 
bool Point3D::operator>=(const Point3D& point) const
{
	return !(*this < point);
}

inline 
bool Point3D::operator!=(const Point3D& point) const
{
	return !(*this == point);
}

inline 
bool Point3D::operator>(const Point3D& point) const
{
	return !(*this <= point);
}

inline
double Point3D::operator*(const Point3D& p) const
{
	return x*p.x + y*p.y + z*p.z;
}

inline
Point3D& Point3D::operator=(const Point3D& rhs)
{
	if (this != &rhs)
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
	}
	return *this;
}

inline
bool Point3D::operator<(const Point3D& point) const
{
	if (TOL_LT(x, point.x, tolx)) 
		return true;
	if (TOL_EQ(x, point.x, tolx) && TOL_LT(y, point.y, toly))
		return true;
	if (TOL_EQ(x, point.x, tolx) && TOL_EQ(y, point.y, toly) && TOL_LT(z, point.z, tolz))
		return true;
	return false;
}

inline
std::ostream& operator<<(std::ostream& theStream, Point3D thePoint)
{
	theStream << "(" << thePoint.x << ", " << thePoint.y << ", " << thePoint.z << ")";
    return theStream;
}

inline
std::istream& operator>>(std::istream& theStream, Point3D& thePoint)
{
  theStream >> thePoint.x >> thePoint.y >> thePoint.z;
    return theStream;
}

inline
Point3D Point3D::crossProduct(const Point3D& p) const
{
	Point3D outCP;
	outCP.x =  y*p.z - z*p.y;
	outCP.y = -x*p.z + z*p.x;
	outCP.z =  x*p.y - y*p.x;
	return outCP;
}

inline
double Point3D::crossProductNorm(const Point3D& point) const
{
	return abs(crossProduct(point));
}

inline
Point3D crossProduct(const Point3D& p1, const Point3D& p2)
{
	Point3D outCP;
	outCP.x =  p1.y*p2.z - p1.z*p2.y;
	outCP.y = -p1.x*p2.z + p1.z*p2.x;
	outCP.z =  p1.x*p2.y - p1.y*p2.x;
	return outCP;
}

inline
double triangleArea(const Point3D& p1, const Point3D& p2, const Point3D& p3)
{
	Point3D v1 = p2 - p1, v2 = p3 - p1;
	Point3D tmp = crossProduct(v1, v2);
	return 0.5*abs(tmp);
}

inline
double tetrahedronVolume(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4)
{
	Point3D rXi1(p2 - p1), rXi2(p3 - p1), rXi3(p4 - p1);
	Point3D temp = rXi2.crossProduct(rXi3);
	return (rXi1*temp)/6.;
}

inline
Point3D centroid(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4)
{
	Point3D outPt;
	outPt.x = 0.25*(p1.x + p2.x + p3.x + p4.x);
	outPt.y = 0.25*(p1.y + p2.y + p3.y + p4.y);
	outPt.z = 0.25*(p1.z + p2.z + p3.z + p4.z);
	return outPt;
}

inline
double solidAngle(const Point3D& o, const Point3D& a, const Point3D& b, const Point3D& c)
{
	Point3D ra(a - o), rb(b - o), rc(c - o);
	Point3D temp = rb.crossProduct(rc);
	double detabc = ra*temp;
	double la = abs(ra), lb = abs(rb), lc = abs(rc);
	double denom = la*lb*lc + lc*(ra*rb) + lb*(ra*rc) + la*(rb*rc);
	double sa = atan(fabs(detabc)/denom);
	if(sa < 0.)
		sa += Topologies::PI;
	return 2.*sa;
}

#endif

