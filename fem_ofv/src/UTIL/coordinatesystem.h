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

#ifndef COORDINATESYSTEM_H
#define COORDINATESYSTEM_H

struct Point3D;

//! Namespace for computing coordinate system conversions
namespace CoordinateSystem{
	enum class Type{cartesian, cylindrical, spherical};

	//! Convert a point in cylindrical coordinates to Cartesian coordinates
	/*! The template parameter specifies either Point2D or Point3D
	*/
	template<typename PointType>
  PointType cylToCart(PointType const& p);
	//! Convert a vector in cylindrical coordinates at a point p (in Cartesian) to a vector in Cartesian coordinates
	/*! The template parameter specifies either Point2D or Point3D
	 *  @param v is a vector in cylindrical coordinates, meaning that
	 *  v.x is v.x \hat{\rho}, v.y is v.y \hat{\phi}, and v.z is just v.z
	*/
	template<typename PointType>
	PointType cylVectorToCartVector(PointType const& v, PointType const& p);
	//! Convert a point in spherical coordinates to Cartesian coordinates
	Point3D sphrToCart(Point3D const& p);
	//! Convert a vector in spherical coordinates at a point p (in Cartesian) to a vector in Cartesian coordinates
	/*! @param v is a vector in spherical coordinates, meaning that 
	 *  v.x is v.x \hat{r}, v.y is v.y \hat{\theta}, and v.z is v.z \hat{\phi}
	*/
	Point3D sphrVectorToCartVector(Point3D const& v, Point3D const& p);
	//! Converts vector @param v in system specified by @param ct to Cartesian coordinates at point @param p
	Point3D convertVector(Point3D const& v, Point3D const& p, Type ct);
}

#endif

