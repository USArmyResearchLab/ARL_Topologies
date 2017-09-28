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

#include "filter3d.h"
#include "tomesh.h"
#include "tomeshprocessing.h"
#include "cartesianmesher.h"
#include <list>
#include <CGAL/Fuzzy_sphere.h>

namespace Topologies{
template <class WeightFunc>
Filter3D<WeightFunc>::Filter3D(const TOMesh* const inMesh, bool ptsAtCentroids) 
{
	assert(inMesh->dimNum() == 3);
	if(ptsAtCentroids)
		setupFromMesh(inMesh);
	else
		setupFromMeshNodal(inMesh);
}

template <class WeightFunc>
Filter3D<WeightFunc>::Filter3D(const TOMesh* const inMesh, double inRad, bool ptsAtCentroids) :
	FilterBase(inRad)
{
	assert(inMesh->dimNum() == 3);
	if(ptsAtCentroids)
		setupFromMesh(inMesh);
	else
		setupFromMeshNodal(inMesh);
}

template <class WeightFunc>
void Filter3D<WeightFunc>::setupFromMesh(const TOMesh* const inMesh)
{
	assert(inMesh->dimNum() == 3);
	std::size_t nelem = inMesh->getNumElements();
	ptVec.reserve(nelem);
	areaVec.reserve(nelem);
	valVec.reserve(nelem);
	for(std::size_t k = 0; k < nelem; ++k)
	{
		ptVec.push_back(TOMeshProcessing::getElementCentroid3D(k, inMesh));
		areaVec.push_back(TOMeshProcessing::computeElementVolume(k, inMesh));
		valVec.push_back(inMesh->getOptVal(k)); // Used for interpolation
	}
	searchTree = std::unique_ptr<CGAL::Kd_tree<Traits>>(new CGAL::Kd_tree<Traits>());
	searchTree->insert(ptVec.begin(), ptVec.end());
	for(std::size_t k = 0; k < ptVec.size(); ++k)
		ptIDMap[ptVec[k]] = k;
}

template <class WeightFunc>
void Filter3D<WeightFunc>::setupFromMeshNodal(const TOMesh* const inMesh)
{
	assert(inMesh->dimNum() == 3);
	std::size_t npts = inMesh->getNumNodes();
	ptVec.reserve(npts);
	areaVec.reserve(npts);
	valVec.reserve(npts);
	for(std::size_t k = 0; k < npts; ++k)
	{
		ptVec.push_back(inMesh->getNode3D(k));
		areaVec.push_back(1.);
		valVec.push_back(1.);
	}
	searchTree = std::unique_ptr<CGAL::Kd_tree<Traits>>(new CGAL::Kd_tree<Traits>());
	searchTree->insert(ptVec.begin(), ptVec.end());
	for(std::size_t k = 0; k < ptVec.size(); ++k)
		ptIDMap[ptVec[k]] = k;
}

template <class WeightFunc>
Filter3D<WeightFunc>::Filter3D(unsigned nx, unsigned ny, unsigned nz, double width, double length, 
	double height, bool ptsAtCentroids)
{
	std::vector<double> optVals(nx*ny*nz,1.);
  std::unique_ptr<TOMesh> tmpMesh = CartesianMesher::generateMesh(metHex, nx, ny, nz, width, length, height, optVals);
	assert(tmpMesh);
	if(ptsAtCentroids)
	  setupFromMesh(tmpMesh.get());
	else
		setupFromMeshNodal(tmpMesh.get());
}

template <class WeightFunc>
Filter3D<WeightFunc>::Filter3D(const Filter3D<WeightFunc>& copy) :
	FilterBase(copy),
	ptVec(copy.ptVec),
	ptIDMap(copy.ptIDMap)
{
	searchTree = std::unique_ptr<CGAL::Kd_tree<Traits>>(new CGAL::Kd_tree<Traits>());
	searchTree->insert(ptVec.begin(), ptVec.end());
}

template <class WeightFunc>
Filter3D<WeightFunc>::Filter3D(Filter3D && copy)
{
	swap(copy);
}

template <class WeightFunc>
Filter3D<WeightFunc>& Filter3D<WeightFunc>::operator=(Filter3D<WeightFunc> copy)
{
	swap(copy);
	return *this;
}

template <class WeightFunc>
void Filter3D<WeightFunc>::swap(Filter3D<WeightFunc>& arg2)
{
	FilterBase::swap(arg2);
	searchTree.swap(arg2.searchTree);
  ptVec.swap(arg2.ptVec);
  ptIDMap.swap(arg2.ptIDMap);
}

template <class WeightFunc>
std::vector<double> Filter3D<WeightFunc>::operator()(const std::vector<double>& xVec, const std::vector<Point_3_base>& filtPts, 
		double rad) const
{
	return filterMesh(xVec, filtPts, rad);
}

template <class WeightFunc>
std::vector<double> Filter3D<WeightFunc>::operator()(const std::vector<double>& xVec, double rad) const
{
	return filterMesh(xVec, ptVec, rad);
}

template <class WeightFunc>
Tr_GT::FT Filter3D<WeightFunc>::operator()(Tr_GT::Point_3 p) const
{
	return filterOnePoint(p, curRad, valVec);
}

template <class WeightFunc>
std::vector<double> Filter3D<WeightFunc>::filterMesh(const std::vector<double>& x, const std::vector<Point_3_base>& filtPts, 
		double rad) const
{
	assert(x.size() == ptVec.size());
	std::vector<double> y(filtPts.size());
	for(std::size_t k = 0; k < filtPts.size(); ++k)
		y[k] = filterOnePoint(filtPts[k], rad, x);
	return y;
}

template <class WeightFunc>
double Filter3D<WeightFunc>::filterOnePoint(const Point_3_base& pt, double rad, const std::vector<double>& x) const
{
	WeightFunc wf(rad);
	std::vector<Point_3_base> ptsInRange;
	getPointsInSphere(pt, rad, ptsInRange);
	if(ptsInRange.empty())
		return 0.;
	double sum = 0.;
	double normVal = 0.;
	for(std::size_t k2 = 0; k2 < ptsInRange.size(); ++k2)
	{
		double d2 = CGAL::squared_distance(pt, ptsInRange[k2]);
		double f = wf(sqrt(d2));
		std::unordered_map<Point_3_base, std::size_t, Point_3_hash>::const_iterator mcit = ptIDMap.find(ptsInRange[k2]);
		if(mcit != ptIDMap.end())
		{
			std::size_t id = mcit->second;
			sum += f*areaVec[id]*x[id];
			normVal += f*areaVec[id];
		}
		else
			std::cout << "Warning: point " << ptsInRange[k2] << " not found in ptIDMap" << std::endl;
	}
	return sum/normVal;
}

template <class WeightFunc>
std::vector<std::map<std::size_t, double>> Filter3D<WeightFunc>::diffFilter(const std::vector<Point_3_base>& filtPts, 
	double rad) const
{
	WeightFunc wf(rad);
	std::vector<std::map<std::size_t, double>> res(ptVec.size());
	for(std::size_t k = 0; k < filtPts.size(); ++k)
	{
		std::vector<Point_3_base> ptsInRange;
		getPointsInSphere(filtPts[k], rad, ptsInRange);
		if(ptsInRange.empty())
			continue;
		// First compute normalization and weight function values
		double normVal = 0.;
		std::vector<double> fVec(ptsInRange.size());
		for(std::size_t k2 = 0; k2 < ptsInRange.size(); ++k2)
		{
			std::unordered_map<Point_3_base, std::size_t, Point_3_hash>::const_iterator mcit = ptIDMap.find(ptsInRange[k2]);
			if(mcit != ptIDMap.end())
			{
				double d2 = CGAL::squared_distance(filtPts[k], ptsInRange[k2]);
				double f = wf(sqrt(d2));
				std::size_t id = mcit->second;
				normVal += f*areaVec[id];
				fVec[k2] = f*areaVec[id];
			}
		}
		// Compute derivative
		for(std::size_t k2 = 0; k2 < ptsInRange.size(); ++k2)
		{
			std::unordered_map<Point_3_base, std::size_t, Point_3_hash>::const_iterator mcit = ptIDMap.find(ptsInRange[k2]);
			if(mcit != ptIDMap.end())
			{
				std::size_t id = mcit->second;
				std::map<std::size_t, double>& curRes = res[id];
				curRes.emplace(k, fVec[k2]/normVal);
			}
			else
				std::cout << "Warning: point " << ptsInRange[k2] << " not found in ptIDMap" << std::endl;
		}
	}
	return res;
}

template <class WeightFunc>
void Filter3D<WeightFunc>::getPointsInSphere(const Point_3_base& center, double rad, std::vector<Point_3_base>& outVec) const
{
	CGAL::Fuzzy_sphere<Traits> fs(center, rad);
	outVec.clear();
	searchTree->search(std::back_inserter(outVec), fs);
}

template class Filter3D<HelperNS::linearHat>;
template class Filter3D<HelperNS::constantFunction>;
}// namespace
