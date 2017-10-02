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

#ifndef FILTER2D_H
#define FILTER2D_H

#include "filterbase.h"
#include "cgal_types.h"
#include "helper.h"
#include <vector>
#include <unordered_map>
#include <map>
#include <CGAL/Point_set_2.h>

namespace Topologies{
class TOMesh;

//! A hash function for Point_2_base
/*! Adapted from http://en.cppreference.com/w/cpp/utility/hash */
struct Point_2_hash
{
	std::size_t operator()(Point_2_base const& s) const
	{
		std::size_t const h1 ( std::hash<double>()(s.x()) );
		std::size_t const h2 ( std::hash<double>()(s.y()) );
		return h1 ^ (h2 << 1);
	}
};

//! A two-dimensional filter for filtering data on a mesh or grid
/*! This class implements a 2d filter with a linear (hat) kernel of a given radius.  It is implemented
 *  using the CGAL class Point_set_2 for fast point searching.
 */
template <class WeightFunc = HelperNS::linearHat>
class Filter2D : public FilterBase
{
public:
	Filter2D();
	//! Constructor taking a TOMesh, which sets up a filter to filter data defined either on the elements or nodes
	/*! Input argument `ptsAtCentroids` specifies either an element-based filter (true) or a nodal based filter (false)
	 */
	Filter2D(const TOMesh* const inMesh, bool ptsAtCentroids = true);
	//! Constructor that sets up a 2D filter on a grid with dimensions nx, ny and physical size width, height
	/*! As with the mesh-based constructor, the data can be defined on either the nodes or element centers */
	Filter2D(unsigned nx, unsigned ny, double width, double height, bool ptsAtCenters = true);
	virtual ~Filter2D() {}

	virtual std::unique_ptr<FilterBase> clone() const {return std::unique_ptr<FilterBase>(new Filter2D<WeightFunc>(*this));}
	virtual std::vector<double> operator()(const std::vector<double>& xVec, double rad) const;
	virtual std::vector<double> operator()(const std::vector<double>& xVec, 
		const std::vector<Point_2_base>& filtPts, double rad) const;
	virtual Tr_GT::FT operator()(Tr_GT::Point_2 p) const;
	virtual std::vector<std::map<std::size_t, double>> diffFilter(const std::vector<Point_2_base>& filtPts, double rad) const;
	
private:
	void filterMesh(const std::vector<double>& x, std::vector<double>& y, const std::vector<Point_2_base>& filtPts, double rad) const;
	void filterMesh(const std::vector<double>& x, std::vector<double>& y, double rad) const;
	void getPointsInCircle(const Point_2_base& center, double rad, std::vector<CGAL::Point_set_2<Mesh_K>::Vertex_handle>& outVec) const;
private:
// Methods
	void finishSetupCentroid(const TOMesh* const inMesh);
	void finishSetupNode(const TOMesh* const inMesh);
	double filterOnePoint(const Point_2_base& pt, double rad, const std::vector<double>& x) const;

// Data
	mutable CGAL::Point_set_2<Mesh_K> pointSet;
	std::vector<Point_2_base> ptVec;
	std::unordered_map<Point_2_base, std::size_t, Point_2_hash> ptIDMap;
	// Pixel based filter
	std::vector<unsigned> sizes;
	std::vector<double> physicalSizes;
private:
	// Functions using 3d points
	virtual std::vector<double> operator()(const std::vector<double>& xVec,
		const std::vector<Point_3_base>& filtPts, double rad) const {return std::vector<double>();}
	virtual Tr_GT::FT operator()(Tr_GT::Point_3 p) const {return 0.;}
	virtual std::vector<std::map<std::size_t, double>> diffFilter(const std::vector<Point_3_base>& filtPts, double rad) const
  {return std::vector<std::map<std::size_t, double>>();}
};
}// namespace
#endif

