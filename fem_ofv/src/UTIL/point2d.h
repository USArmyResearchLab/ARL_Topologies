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

#ifndef POINT2D_H
#define POINT2D_H

#include "UTIL/topologiesdefs.h"
#include "REP/cgal_types.h"
#include <iostream>
#include <assert.h>

//! A basic implementation of a point in 2d space
/*! The class has publically accesible coordinates (x,y) and numerous functions
 *  to perform operations on 2d points.  
 */
struct Point2D
{
public:
	//! Point coordinates
	double x, y;
	//! Default constructor initializes all coordinates to 0
	Point2D();
	//! Constructor that sets all coordinates to the given arguments
	Point2D(double inX, double inY);
	//! Constructor that sets all coordinates to the given arguments, and ignores bogusZ
	/*! This is useful for template functions w/ Point2D and Point3D as templates*/
	Point2D(double inX, double inY, double bogusZ);
	Point2D(const Point2D& point); 
	Point2D(const Point2D* point);
	//! Converts a CGAL Point_2 to a Point2D
	Point2D(const Topologies::Point_2_base& point);
	Point2D& operator=(const Point2D& rhs);
	//! Returns the negation of this point, i.e. (-x, -y)
	Point2D operator-() const;
	//! Returns the sum of this point and the argument p
	Point2D operator+(const Point2D& p) const;
	//! Returns the difference of this point and the argument p
	Point2D operator-(const Point2D& p) const;
	//! Sums the argument and stores in this Point3D
	const Point2D& operator+=(const Point2D& p);
	//! Computes the difference of the argument and stores result in this Point3D
	const Point2D& operator-=(const Point2D& p);
	//! Scalar multiplication of a Point2D
	friend Point2D operator*(double a, const Point2D& p);
	//! Scalar multiplication of a Point2D
	friend Point2D operator*(const Point2D& p, double a); 
	//! Scalar product (dot product) of two Point2D objects
	double operator*(const Point2D& p) const;
	//! Computes and stores the result of scalar multiplication with argument a
	const Point2D& operator*=(double a);
	//! Scalar division
	Point2D operator/(double a) const;
	//! Computes and stores the result of scalar division with argument a
	const Point2D& operator/=(double a);
	//! Returns the length of the vector connecting argument v and the origin
	/*! Gives \rho in polar coordinates */
	friend double abs(const Point2D& v);
	//! Returns the angle of the vector connecting argument v and the origin with respect to the x-axis
	/*! Gives the angle in polar coordinates*/
	friend double arg(const Point2D& v);
	//! Checks lexographical order of this Point2D and the argument.
	/*! Ordering is given by first checking x-coordinate.  If the x-coordinates are equal, the y-coordinate is used.*/
	bool operator<(const Point2D& point) const;
	//! Checks equality of the point's coordinates
	bool operator==(const Point2D& point) const;
	//! Checks lexographical order of this Point2D and the argument.
	bool operator<=(const Point2D& point) const;
	//! Checks lexographical order of this Point2D and the argument.
	bool operator>=(const Point2D& point) const;
	//! Checks equality of the point's coordinates
	bool operator!=(const Point2D& point) const;
	//! Checks lexographical order of this Point2D and the argument.
	bool operator>(const Point2D& point) const;
	//! Checks if the arguments are in counter-clockwise order based on their angles
	bool ccwsort(const Point2D& p1, const Point2D& p2) const;
	//! Text output to the specified output stream.  Outputs the coordinate values
	friend std::ostream& operator<<(std::ostream& theStream, Point2D thePoint);
	//! Text input from the specified input stream.  Sets all three coordinates
	friend std::istream& operator>>(std::istream& theStream, Point2D& thePoint);
	//! Computes and returns the vector norm of the cross product of this Point2D and the argument
	double crossProductNorm(const Point2D& point) const;
	//! Returns the dimension of the space this object lives in
	unsigned dim() const {return 2;}
	//! Coordinate access (const)
	double operator[](unsigned k) const;
private:
	static double tol, tolx, toly;
};

inline
double Point2D::operator[](unsigned k) const
{
	if(k == 0) 
		return x;
	else if(k == 1)
		return y;
	return 0.;
}


inline
Point2D::Point2D():
	x(0.), y(0.)
{
}

inline
Point2D::Point2D(double inX, double inY):
	x(inX), y(inY)
{
}

inline
Point2D::Point2D(const Topologies::Point_2_base& point) :
	x(point.x()), y(point.y())
{
}

inline
Point2D::Point2D(double inX, double inY, double bogusZ):
	x(inX), y(inY)
{
}

inline 
Point2D::Point2D(const Point2D& point):
	x(point.x), y(point.y)
{
}
 
inline 
Point2D::Point2D(const Point2D* point):
x(point->x), y(point->y)
{
}

inline
Point2D Point2D::operator-() const
{
	return Point2D(-x, -y);
}

inline
Point2D Point2D::operator+(const Point2D& p) const
{
	return Point2D(x + p.x, y + p.y);
}

inline
Point2D Point2D::operator-(const Point2D& p) const
{
	return Point2D(x - p.x, y - p.y);
}

inline 
const Point2D& Point2D::operator+=(const Point2D& p)
{
	x += p.x;
	y += p.y;
	return *this;
}

inline 
const Point2D& Point2D::operator-=(const Point2D& p)
{
	x -= p.x;
	y -= p.y;
	return *this;
}

inline
Point2D operator*(double a, const Point2D& p)
{
	return Point2D(a*p.x, a*p.y);
}

inline
Point2D operator*(const Point2D& p, double a)
{
	return Point2D(a*p.x, a*p.y);
}

inline 
const Point2D& Point2D::operator*=(double a)
{
	x *= a;
	y *= a;
	return *this;
}

inline
Point2D Point2D::operator/(double a) const
{
	return Point2D(x/a, y/a);
}

inline 
const Point2D& Point2D::operator/=(double a)
{
	x /= a;
	y /= a;
	return *this;
}

inline 
double abs(const Point2D& v)
{
	return std::sqrt(v.x*v.x + v.y*v.y);
}

inline 
double arg(const Point2D& v)
{
	return atan2(v.y, v.x);
}

inline 
bool Point2D::operator==(const Point2D& point) const
{
	return (TOL_EQ(x, point.x, tolx)) && (TOL_EQ(y, point.y, toly));
}

inline 
bool Point2D::operator<=(const Point2D& point) const
{
	return (*this < point) || (*this == point);
}

inline 
bool Point2D::operator>=(const Point2D& point) const
{
	return !(*this < point);
}

inline 
bool Point2D::operator!=(const Point2D& point) const
{
	return !(*this == point);
}

inline 
bool Point2D::operator>(const Point2D& point) const
{
	return !(*this <= point);
}

inline
double Point2D::operator*(const Point2D& p) const
{
	return x*p.x + y*p.y;
}

inline
double Point2D::crossProductNorm(const Point2D& point) const
{
	return x*point.y - y*point.x;
}

#endif

