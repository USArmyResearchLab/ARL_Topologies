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

#ifndef FILTERBASE_H
#define FILTERBASE_H

#include "cgal_types.h"
#include "helper.h"
#include <vector>
#include <memory>

namespace Topologies{
//! Abstract base class for 2d or 3d filters
class FilterBase
{
	public:
	public:
	FilterBase() : curRad(0.1) {}
	explicit FilterBase(double inRad) : curRad(inRad) {}
	virtual ~FilterBase() {}
	//! Clone function, ie. polymorphic copy
	virtual std::unique_ptr<FilterBase> clone() const = 0;
	//! Returns a vector of filtered data using a filter of radius `rad`, where xVec is assumed to lie on the previously specified grid or mesh
	virtual std::vector<double> operator()(const std::vector<double>& xVec, double rad) const = 0;
	//! Returns a vector of filtered data using a filter of radius `rad`, where xVec is assumed to coincide with the points specified in the original mesh, and the results are given on filtPts
	virtual std::vector<double> operator()(const std::vector<double>& xVec, 
		const std::vector<Point_2_base>& filtPts, double rad) const = 0;
	//! Returns a vector of filtered data using a filter of radius `rad`, where xVec is assumed to coincide with the points specified in the original mesh, and the results are given on filtPts
	virtual std::vector<double> operator()(const std::vector<double>& xVec,
		const std::vector<Point_3_base>& filtPts, double rad) const = 0;
	//! Returns a sparse-matrix like vector of vectors of column/value pairs that gives the partial derivatives of values at filtPts, with respect to all of the original nodes used to initialize the FilterBase object.
	/*! At return, the vector has size of the number of original points (nodes in the input mesh or grid)
	 *  and the vector associated with that index contains all non-zero partial derivatives at parameter filtPts. 
	 */
	virtual HelperNS::SparseMatrix diffFilter(const std::vector<Point_2_base>& filtPts, double rad) const = 0;
	//! Returns a sparse-matrix like vector of vectors of column/value pairs that gives the partial derivatives of values at filtPts, with respect to all of the original nodes used to initialize the FilterBase object.
	virtual HelperNS::SparseMatrix diffFilter(const std::vector<Point_3_base>& filtPts, double rad) const = 0;

	//! Returns the filtered value of the data loaded in the TOMesh at object construction, at the point `p`
	virtual Tr_GT::FT operator()(Tr_GT::Point_2 p) const = 0;
	//! Returns the filtered value of the data loaded in the TOMesh at object construction, at the point `p`
	virtual Tr_GT::FT operator()(Tr_GT::Point_3 p) const = 0;
	//! Set the current filter radius
	void setCurRad(double inRad) {curRad = inRad;}
	//! Set the vector of values used for filtering
	void setValVec(const std::vector<double>& inVV) {valVec = inVV;}
protected:
	void swap(FilterBase& arg2);

	double curRad;
	std::vector<double> areaVec, valVec;
};

inline
void FilterBase::swap(FilterBase& arg2)
{
	std::swap(curRad, arg2.curRad);
	areaVec.swap(arg2.areaVec);
	valVec.swap(arg2.valVec);
}
}//namespace
#endif

