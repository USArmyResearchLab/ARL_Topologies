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

// Questions?
// Contact: Raymond Wildman, raymond.a.wildman.civ@mail.mil

#ifndef FILTER3D_H
#define FILTER3D_H

#include "filterbase.h"
#include "cgal_types.h"
#include "helper.h"
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>

namespace Topologies{
class TOMesh;

//! A hash function for Point_3_base
/*! Adapted from http://en.cppreference.com/w/cpp/utility/hash */
struct Point_3_hash
{
	std::size_t operator()(Point_3_base const& s) const
	{
		std::size_t const h1 ( std::hash<double>()(s.x()) );
		std::size_t const h2 ( std::hash<double>()(s.y()) );
		std::size_t const h3 ( std::hash<double>()(s.z()) );
		return h1 ^ ((h2 ^ (h3 << 1)) << 1);
	}
};

//! A three-dimensional filter for filtering data on a mesh or grid
/*! This class implements a 3d filter with a linear (hat) kernel of a given radius.  It is implemented
 *  using the CGAL class Kd_tree for fast point searching.
 */
template <class WeightFunc = HelperNS::linearHat>
class Filter3D : public FilterBase
{
public:
	//! Constructor taking a TOMesh as an argument, sets up a filter assuming data will be located at element centroids
	/*! Data contained in the TOMesh as `optVal` will also be loaded and can be filtered using operator() */
	Filter3D(const TOMesh* const inMesh, bool ptsAtCentroids = true);
	//! Constructor taking a TOMesh as an argument and a filter radius, sets up a filter assuming data will be located at element centroids
	/*! Data contained in the TOMesh as `optVal` will also be loaded and can be filtered using operator() */
	Filter3D(const TOMesh* const inMesh, double inRad, bool ptsAtCentroids = true);
	//! Constructor to set up a filter on a structured grid with dimensions nx, ny, nz, and size width, length, height
	Filter3D(unsigned nx, unsigned ny, unsigned nz, double width, double length, double height, bool ptsAtCentroids = true);
	virtual ~Filter3D() {}

	Filter3D(const Filter3D& copy);
	Filter3D(Filter3D && copy);
	Filter3D& operator=(Filter3D copy);
	void swap(Filter3D& arg2);
	
	virtual std::unique_ptr<FilterBase> clone() const {return std::unique_ptr<FilterBase>(new Filter3D<WeightFunc>(*this));}
	virtual std::vector<double> operator()(const std::vector<double>& xVec, double rad) const;
	virtual std::vector<double> operator()(const std::vector<double>& xVec, 
		const std::vector<Point_3_base>& filtPts, double rad) const;
	virtual Tr_GT::FT operator()(Tr_GT::Point_3 p) const;
	virtual std::vector<std::map<std::size_t, double>> diffFilter(const std::vector<Point_3_base>& filtPts, double rad) const;
private:
	std::vector<double> filterMesh(const std::vector<double>& x, const std::vector<Point_3_base>& filtPts, double rad) const;
	void getPointsInSphere(const Point_3_base& center, double rad, std::vector<Point_3_base>& outVec) const;
	double filterOnePoint(const Point_3_base& pt, double rad, const std::vector<double>& y) const;
	void setupFromMesh(const TOMesh* const inMesh);
	void setupFromMeshNodal(const TOMesh* const inMesh);
	double weightFunc(double d, double rad) const;
private:
	typedef CGAL::Search_traits_3<Mesh_K> Traits;
	std::unique_ptr<CGAL::Kd_tree<Traits>> searchTree;
	std::vector<Point_3_base> ptVec;
	std::unordered_map<Point_3_base, std::size_t, Point_3_hash> ptIDMap;
private:
	virtual Tr_GT::FT operator()(Tr_GT::Point_2 p) const {return 0.;}
	virtual std::vector<double> operator()(const std::vector<double>& xVec, 
		const std::vector<Point_2_base>& filtPts, double rad) const {return std::vector<double>();}
	virtual std::vector<std::map<std::size_t, double>> diffFilter(const std::vector<Point_2_base>& filtPts, double rad) const
	{return std::vector<std::map<std::size_t, double>>();}
};
} //namespace
#endif
