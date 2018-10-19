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

#define CATCH_CONFIG_MAIN

#include "filter2d.h"
#include "tomesh.h"
#include "tomeshprocessing.h"
#include "catch.hpp"
#include <iostream>
#include <fstream>

using namespace Topologies;
using namespace HelperNS;

std::vector<double> filterWithDiff(const std::vector<double>& xVec, std::size_t npts,
		SparseMatrix const& diffMat)
{
	REQUIRE(xVec.size() == diffMat.size());
	std::vector<double> res = diffMat.transposeTimes(xVec, npts);
	return res;
}

TEST_CASE("Testing Filter2D with uniform mesh input.", "[Filter2D]")
{
	unsigned nx = 20, ny = 20;
	double width = 1., height = 1.;
	bool ptsAtCenter = true;
	SECTION("Testing centroid points")
	{
		Filter2D<> testFilt(nx, ny, width, height, ptsAtCenter);
		std::vector<double> xVec(nx*ny);
		double rad = 0.2;
		for(unsigned ky = 0; ky < ny; ++ky)
		{
			for(unsigned kx = 0; kx < nx; ++kx)
			{
				if(kx < nx/2)
					xVec[ky*nx + kx] = 1.;
				else
					xVec[ky*nx + kx] = 0.;
			}
		}
		std::vector<double> res = testFilt(xVec, rad);
		REQUIRE(res.size() == xVec.size());
		// Check middle row
		std::vector<double> row = {1., 1., 1., 1., 1., 1., 1., 0.9482920597029305, 0.8190363597341801, 0.6194060825994609, 
			0.3805939174005385, 0.180963640265819, 0.05170794029706944, 0., 0., 0., 0., 0., 0., 0.};
		for(unsigned kx = 0; kx < nx; ++kx)
			REQUIRE(res[(ny/2)*nx + kx] == Approx(row[kx]));
		// Check symmetry
		for(unsigned ky = 0; ky < ny; ++ky)
			for(unsigned kx = 0; kx < nx/2; ++kx)
				REQUIRE(fabs(res[ky*nx + kx] - 0.5) == Approx(fabs(res[ky*nx + nx - kx - 1] - 0.5)));
		for(unsigned ky = 0; ky < ny/2; ++ky)
			for(unsigned kx = 0; kx < nx; ++kx)
				REQUIRE(res[ky*nx + kx] == Approx(res[(ny - ky - 1)*nx + kx]));
		// Test single points
		std::vector<Point_2_base> ptVec = {Point_2_base(0., 0.5), Point_2_base(0.5, 0.5), Point_2_base(1., 0.5), Point_2_base(2., 2.)};
		res = testFilt(xVec, ptVec, rad);
		REQUIRE(res.size() == ptVec.size()); // Check sizes
		row = {1., 0.5, 0., 0.}; // Check values
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
		// Retest with the filter derivative (same results can be computed since it's a linear filter)
		SparseMatrix diffMat = testFilt.diffFilter(ptVec, rad);
		res = filterWithDiff(xVec, ptVec.size(), diffMat);
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
	}
	ptsAtCenter = false;
	SECTION("Testing points at nodes")
	{
		Filter2D<> testFilt0(nx, ny, width, height, ptsAtCenter);
		// use copy ctor
		Filter2D<> testFilt(testFilt0);
		double rad = 0.2;
		++nx;
		++ny;
		std::vector<double> xVec(nx*ny);
		for(unsigned ky = 0; ky < ny; ++ky)
		{
			for(unsigned kx = 0; kx < nx; ++kx)
			{
				if(kx < nx/2)
					xVec[ky*nx + kx] = 1.;
				else if(kx == nx/2)
					xVec[ky*nx + kx] = 0.5;
				else
					xVec[ky*nx + kx] = 0.;
			}
		}
		std::vector<double> res = testFilt(xVec, rad);
		REQUIRE(res.size() == xVec.size());
		// Check middle row
		std::vector<double> row = {1., 1., 1., 1., 1., 1., 1., 0.9741460298514653, 0.8836642097185555, 0.7192212211668211, 
			0.500000000000000, 0.2807787788331791, 0.1163357902814442, 0.0258539701485347, 0., 0., 0., 0., 0., 0., 0.};
		std::cout.precision(16);
		for(unsigned kx = 0; kx < nx; ++kx)
			REQUIRE(res[(ny/2)*nx + kx] == Approx(row[kx]));
		// Check symmetry
		for(unsigned ky = 0; ky < ny; ++ky)
			for(unsigned kx = 0; kx < nx/2; ++kx)
				REQUIRE(fabs(res[ky*nx + kx] - 0.5) == Approx(fabs(res[ky*nx + nx - kx - 1] - 0.5)));
		for(unsigned ky = 0; ky < ny/2; ++ky)
			for(unsigned kx = 0; kx < nx; ++kx)
				REQUIRE(res[ky*nx + kx] == Approx(res[(ny - ky - 1)*nx + kx]));
		std::vector<Point_2_base> ptVec = {Point_2_base(0., 0.5), Point_2_base(0.5, 0.5), Point_2_base(1., 0.5), Point_2_base(2., 2.)};
		res = testFilt(xVec, ptVec, rad);
		REQUIRE(res.size() == ptVec.size()); // Check sizes
		row = {1., 0.5, 0., 0.}; // Check values
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
		// Retest with the filter derivative (same results can be computed since it's a linear filter)
		SparseMatrix diffMat = testFilt.diffFilter(ptVec, rad);
		res = filterWithDiff(xVec, ptVec.size(), diffMat);
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
	}
}

TEST_CASE("Testing Filter2D with TOMesh tri input","[Filter2D]")
{
	// Read mesh from file
	std::ifstream meshFile("filter2dmesh.txt");
	REQUIRE(meshFile.good());
	// Read points
	std::size_t npts;
	meshFile >> npts;
	std::vector<Point_2_base> ptVec(npts);
	for(auto pit = ptVec.begin(); pit != ptVec.end(); ++pit)
	{
		double x, y;
		meshFile >> x >> y;
		*pit = Point_2_base(x, y);
	}
	// Read element connectivity
	std::size_t nelems;
	meshFile >> nelems;
	std::vector<std::vector<std::size_t>> connectivity(nelems);
	for(auto cit = connectivity.begin(); cit != connectivity.end(); ++cit)
	{
		std::size_t n1, n2, n3;
		meshFile >> n1 >> n2 >> n3;
		*cit = {n1, n2, n3};
	}
	// Read previously generated results
	std::vector<double> centroidRes(nelems);
	std::size_t ncentres;
	meshFile >> ncentres;
	REQUIRE(ncentres == nelems);
	for(auto cit = centroidRes.begin(); cit != centroidRes.end(); ++cit)
		meshFile >> *cit;
	std::vector<double> nodeRes(nelems);
	std::size_t nnoderes;
	meshFile >> nnoderes;
	REQUIRE(nnoderes == npts);
	for(auto nit = nodeRes.begin(); nit != nodeRes.end(); ++nit)
		meshFile >> *nit;
	// Form mesh and test
	TOMesh2D testMesh(ptVec, connectivity);
	bool ptsAtCenter = true;
	SECTION("Testing centroid points")
	{
		Filter2D<> testFilt0(&testMesh, ptsAtCenter);
		// Use assignment operator
		Filter2D<> testFilt = testFilt0;
		std::vector<double> xVec(testMesh.getNumElements());
		for(std::size_t k = 0; k < xVec.size(); ++k)
		{
			Point_2_base centroid = TOMeshProcessing::getElementCentroid2D(k, &testMesh);
			if(centroid.y() < 0.5)
				xVec[k] = 0.;
			else
				xVec[k] = 1.;
		}
		double rad = 0.15;
		std::vector<double> res = testFilt(xVec, rad);
		REQUIRE(res.size() == xVec.size());
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(centroidRes[k]));
		// Test outside region
		std::vector<Point_2_base> ptVec2 = {Point_2_base(-1., -1.), Point_2_base(-1., 1.), Point_2_base(2., 2.)};
		res = testFilt(xVec, ptVec2, rad);
		REQUIRE(res.size() == ptVec2.size()); // Check sizes
		std::vector<double> row(res.size(), 0.); // Check values
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
		// Retest with the filter derivative (same results can be computed since it's a linear filter)
		SparseMatrix diffMat = testFilt.diffFilter(ptVec2, rad);
		res = filterWithDiff(xVec, ptVec2.size(), diffMat);
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
	}
	ptsAtCenter = false;
	SECTION("Testing points at nodes")
	{
		Filter2D<> testFilt(&testMesh, ptsAtCenter);
		std::vector<double> xVec(testMesh.getNumNodes());
		for(std::size_t k = 0; k < xVec.size(); ++k)
		{
			Point_2_base curNode = testMesh.getNode2D(k);
			if(curNode.y() < 0.5)
				xVec[k] = 0.;
			else
				xVec[k] = 1.;
		}
		double rad = 0.15;
		std::vector<double> res = testFilt(xVec, rad);
		REQUIRE(res.size() == xVec.size());
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(nodeRes[k]));
		// Test outside region
		std::vector<Point_2_base> ptVec2 = {Point_2_base(-1., -1.), Point_2_base(-1., 1.), Point_2_base(2., 2.)};
		res = testFilt(xVec, ptVec2, rad);
		REQUIRE(res.size() == ptVec2.size()); // Check sizes
		std::vector<double> row(res.size(), 0.); // Check values
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
		// Retest with the filter derivative (same results can be computed since it's a linear filter)
		SparseMatrix diffMat = testFilt.diffFilter(ptVec2, rad);
		res = filterWithDiff(xVec, ptVec2.size(), diffMat);
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
	}
}

