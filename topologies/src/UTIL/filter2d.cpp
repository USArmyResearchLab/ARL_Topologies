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

#include "filter2d.h"
#include "helper.h"
#include "tomesh.h"
#include "tomeshprocessing.h"
#include "cartesianmesher.h"
#include <list>

namespace Topologies{
typedef CGAL::Point_set_2<Mesh_K>::Vertex_handle PS_Vertex_handle;

template <class WeightFunc>
Filter2D<WeightFunc>::Filter2D(const TOMesh* const inMesh, bool ptsAtCentroids)
{
	assert(inMesh->dimNum() == 2);
	// Set up CGAL Point_set_2, allows for fast nearest neighbor search
	if(ptsAtCentroids)
		finishSetupCentroid(inMesh);
	else
		finishSetupNode(inMesh);
}

template <class WeightFunc>
Filter2D<WeightFunc>::Filter2D(unsigned nx, unsigned ny, double width, double height, bool ptsAtCentroids)
{
	std::vector<double> optVals(nx*ny,1.);
	std::unique_ptr<TOMesh> tmpMesh = CartesianMesher::generateMesh(metQuad, nx, ny, width, height, optVals);
	if(ptsAtCentroids)
		finishSetupCentroid(tmpMesh.get());
	else
		finishSetupNode(tmpMesh.get());
}

template <class WeightFunc>
void Filter2D<WeightFunc>::finishSetupCentroid(const TOMesh* const inMesh)
{
	// Unknowns associated with elements, use centroids of mesh elements
	std::size_t nelem = inMesh->getNumElements();
	ptVec.reserve(nelem);
	areaVec.reserve(nelem);
	valVec.reserve(nelem);
	for(std::size_t k = 0; k < nelem; ++k)
	{
		ptVec.push_back(TOMeshProcessing::getElementCentroid2D(k, inMesh));
		areaVec.push_back(TOMeshProcessing::computeElementVolume(k, inMesh));
		valVec.push_back(inMesh->getOptVal(k));
	}
	// Insert points into CGAL::Point_set_2 for fast queries
	pointSet.insert(ptVec.begin(), ptVec.end());
	// Set up an unordered map to associate points to ids
	// Necessary since we'll only get a set of points back from Point_set_2
	// but we'll need their indexes into a vector
	for(std::size_t k = 0; k < ptVec.size(); ++k)
		ptIDMap[ptVec[k]] = k;
}

template <class WeightFunc>
void Filter2D<WeightFunc>::finishSetupNode(const TOMesh* const inMesh)
{
	// Unknowns associated with mesh nodes
	std::size_t npts = inMesh->getNumNodes();
	ptVec.reserve(npts);
	areaVec.reserve(npts);
	valVec.reserve(npts);
	for(std::size_t k = 0; k < npts; ++k)
	{
		ptVec.push_back(inMesh->getNode2D(k));
		areaVec.push_back(1.); // Using 1 for now, may be updated for more accurate results
		valVec.push_back(1.);
	}
	// Insert points into CGAL::Point_set_2 for fast queries
	pointSet.insert(ptVec.begin(), ptVec.end());
	// Set up an unordered map to associate points to ids
	// Necessary since we'll only get a set of points back from Point_set_2
	// but we'll need their indexes into a vector
	for(std::size_t k = 0; k < ptVec.size(); ++k)
		ptIDMap[ptVec[k]] = k;
}

template <class WeightFunc>
std::vector<double> Filter2D<WeightFunc>::operator()(const std::vector<double>& xVec, double rad) const
{
	std::vector<double> y(xVec.size());
	filterMesh(xVec, y, rad);
  return y;
}

template <class WeightFunc>
std::vector<double> Filter2D<WeightFunc>::operator()(const std::vector<double>& xVec, const std::vector<Point_2_base>& filtPts, double rad) const
{
	std::vector<double> y(filtPts.size());
	filterMesh(xVec, y, filtPts, rad);
	return y;
}

template <class WeightFunc>
void Filter2D<WeightFunc>::filterMesh(const std::vector<double>& x, std::vector<double>& y, double rad) const
{
	filterMesh(x, y, ptVec, rad);
}

template <class WeightFunc>
void Filter2D<WeightFunc>::filterMesh(const std::vector<double>& x, std::vector<double>& y, const std::vector<Point_2_base>& filtPts, double rad) const
{
	assert(x.size() == ptVec.size());
	if(y.size() != filtPts.size())
		y.resize(filtPts.size());
	for(std::size_t k = 0; k < filtPts.size(); ++k)
		y[k] = filterOnePoint(filtPts[k], rad, x);
}

template <class WeightFunc>
double Filter2D<WeightFunc>::filterOnePoint(const Point_2_base& pt, double rad, const std::vector<double>& x) const
{
	WeightFunc wf(rad);
	std::vector<PS_Vertex_handle> ptsInRange;
	getPointsInCircle(pt, rad, ptsInRange);
	double sum = 0.;
	double normVal = 0.;
	for(std::size_t k2 = 0; k2 < ptsInRange.size(); ++k2)
	{
		double d2 = CGAL::squared_distance(pt, ptsInRange[k2]->point());
		double f = wf(sqrt(d2));
		std::unordered_map<Point_2_base, std::size_t, Point_2_hash>::const_iterator mcit = ptIDMap.find(ptsInRange[k2]->point());
		if(mcit != ptIDMap.end())
		{
			std::size_t id = mcit->second;
			sum += f*areaVec[id]*x[id];
			normVal += f*areaVec[id];
		}
		else
			std::cout << "Warning: point " << ptsInRange[k2]->point() << " not found in ptIDMap" << std::endl;
	}
	if(ptsInRange.empty())
		return 0.;
	return sum/normVal;
}

template <class WeightFunc>
Tr_GT::FT Filter2D<WeightFunc>::operator()(Tr_GT::Point_2 p) const
{
	return filterOnePoint(p, curRad, valVec);
}

template <class WeightFunc>
std::vector<std::map<std::size_t, double>> Filter2D<WeightFunc>::diffFilter(const std::vector<Point_2_base>& filtPts, double rad) const
{
	WeightFunc wf(rad);
	std::vector<std::map<std::size_t, double>> res(ptVec.size());
	for(std::size_t k = 0; k < filtPts.size(); ++k)
	{
		std::vector<PS_Vertex_handle> ptsInRange;
		getPointsInCircle(filtPts[k], rad, ptsInRange);
		// First compute normalization
		double normVal = 0.;
		std::vector<double> fVec(ptsInRange.size());
		for(std::size_t k2 = 0; k2 < ptsInRange.size(); ++k2)
		{
			std::unordered_map<Point_2_base, std::size_t, Point_2_hash>::const_iterator mcit = ptIDMap.find(ptsInRange[k2]->point());
			if(mcit != ptIDMap.end())
			{
				double d2 = CGAL::squared_distance(filtPts[k], ptsInRange[k2]->point());
				double f = wf(sqrt(d2));
				std::size_t id = mcit->second;
				normVal += f*areaVec[id];
				fVec[k2] = f*areaVec[id];
			}
		}
		// Compute derivative
		for(std::size_t k2 = 0; k2 < ptsInRange.size(); ++k2)
		{
			std::unordered_map<Point_2_base, std::size_t, Point_2_hash>::const_iterator mcit = ptIDMap.find(ptsInRange[k2]->point());
			if(mcit != ptIDMap.end())
			{
				std::size_t id = mcit->second;
				std::map<std::size_t, double>& curRes = res[id];
				curRes.emplace(k, fVec[k2]/normVal);
			}
			else
				std::cout << "Warning: point " << ptsInRange[k2]->point() << " not found in ptIDMap" << std::endl;
		}
	}
	return res;
}

template <class WeightFunc>
void Filter2D<WeightFunc>::getPointsInCircle(const Point_2_base& center, double rad, std::vector<PS_Vertex_handle>& outVec) const
{
	CGAL::Circle_2<Mesh_K> rc(center, rad*rad);
	outVec.clear();
	pointSet.range_search(rc, std::back_inserter(outVec));
}

template class Filter2D<HelperNS::linearHat>;
template class Filter2D<HelperNS::constantFunction>;
}//namespace
