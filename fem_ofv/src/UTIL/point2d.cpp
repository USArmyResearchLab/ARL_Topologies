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

#include "point2d.h"

double Point2D::tol = 1e-14, Point2D::tolx = 1e-14, Point2D::toly = 1e-14;

Point2D& Point2D::operator=(const Point2D& rhs)
{
	if (this != &rhs)
	{
		x = rhs.x;
		y = rhs.y;
	}
	return *this;
}

bool Point2D::operator<(const Point2D& point) const
{
	if (TOL_LT(x, point.x, tolx)) 
		return true;
	if (TOL_EQ(x, point.x, tolx) && TOL_LT(y, point.y, toly))
		return true;
	return false;
}

std::ostream& operator<<(std::ostream& theStream, Point2D thePoint)
{
	theStream << "(" << thePoint.x << ", " << thePoint.y << ")";
    return theStream;
}

std::istream& operator>>(std::istream& theStream, Point2D& thePoint)
{
	theStream >> thePoint.x >> thePoint.y; 
    return theStream;
}
