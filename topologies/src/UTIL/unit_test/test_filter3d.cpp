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

#include "filter3d.h"
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

TEST_CASE("Testing Filter3D with uniform mesh input.", "[Filter3D]")
{
	unsigned nx = 20, ny = 20, nz = 20;
	double width = 1., length = 1., height = 1.;
	SECTION("Testing centroid points")
	{
		std::vector<double> xVec(nx*ny*nz);
		for(unsigned kz = 0; kz < nz; ++kz)
		{
			for(unsigned ky = 0; ky < ny; ++ky)
			{
				for(unsigned kx = 0; kx < nx; ++kx)
				{
					if(kx < nx/2)
						xVec[kz*ny*nz + ky*nx + kx] = 1.;
					else
						xVec[kz*ny*nz + ky*nx + kx] = 0.;
				}
			}
		}
		Filter3D<> testFilt(nx, ny, nz, width, length, height, true);
		double rad = 0.2;
		testFilt.setCurRad(rad);
		testFilt.setValVec(xVec);
		std::vector<double> res = testFilt(xVec, rad);
		REQUIRE(res.size() == xVec.size());
		// Check middle row
		std::vector<double> row = {1., 1., 1., 1., 1., 1., 1., 0.960180539842167, 0.8349864995396327, 0.6243290624690414, 
			0.3756709375309588, 0.1650135004603672, 0.03981946015783239, 0., 0., 0., 0., 0., 0., 0.};
		for(unsigned kx = 0; kx < nx; ++kx)
			REQUIRE(res[(nz/2)*ny*nz + (ny/2)*nx + kx] == Approx(row[kx]));
		// Check symmetry
		for(unsigned kz = 0; kz < nz; ++kz)
			for(unsigned ky = 0; ky < ny; ++ky)
				for(unsigned kx = 0; kx < nx/2; ++kx)
					REQUIRE(fabs(res[kz*ny*nz + ky*nx + kx] - 0.5) == Approx(fabs(res[kz*ny*nz + ky*nx + nx - kx - 1] - 0.5)));
		for(unsigned kz = 0; kz < nz; ++kz)
			for(unsigned ky = 0; ky < ny/2; ++ky)
				for(unsigned kx = 0; kx < nx; ++kx)
					REQUIRE(res[kz*ny*nz + ky*nx + kx] == Approx(res[kz*ny*nz + (ny - ky - 1)*nx + kx]));
		for(unsigned kz = 0; kz < nz/2; ++kz)
			for(unsigned ky = 0; ky < ny; ++ky)
				for(unsigned kx = 0; kx < nx; ++kx)
					REQUIRE(res[kz*ny*nz + ky*nx + kx] == Approx(res[(nz - kz - 1)*ny*nz + ky*nx + kx]));
		// Test single points
		std::vector<Point_3_base> ptVec = {Point_3_base(0., 0.5, 0.5), Point_3_base(0.5, 0.5, 0.5), 
			Point_3_base(1., 0.5, 0.5), Point_3_base(2., 2., 2.)};
		row = {1., 0.5, 0., 0.}; // Check values
		for(std::size_t k = 0; k < ptVec.size(); ++k)
			REQUIRE(testFilt(ptVec[k]) == Approx(row[k]));
		// Retest with the filter derivative (same results can be computed since it's a linear filter)
		SparseMatrix diffMat = testFilt.diffFilter(ptVec, rad);
		res = filterWithDiff(xVec, ptVec.size(), diffMat);
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
	}
	SECTION("Testing nodal points")
	{
		Filter3D<> testFilt0(nx, ny, nz, width, length, height, false);
		// Use copy ctor
		Filter3D<> testFilt(testFilt0);
		++nx;
		++ny;
		++nz;
		std::vector<double> xVec(nx*ny*nz);
		for(unsigned kz = 0; kz < nz; ++kz)
		{
			for(unsigned ky = 0; ky < ny; ++ky)
			{
				for(unsigned kx = 0; kx < nx; ++kx)
				{
					if(kx < nx/2)
						xVec[kz*ny*nz + ky*nx + kx] = 1.;
					else if(kx == nx/2)
						xVec[kz*ny*nz + ky*nx + kx] = 0.5;
					else
						xVec[kz*ny*nz + ky*nx + kx] = 0.;
				}
			}
		}
		double rad = 0.2;
		testFilt.setCurRad(rad);
		testFilt.setValVec(xVec);
		std::vector<double> res = testFilt(xVec, rad);
		REQUIRE(res.size() == xVec.size());
		// Check middle row
		std::vector<double> row = {1., 1., 1., 1., 1., 1., 1., 0.9800902699210855, 0.8975835196909002, 0.7296577810043378, 
			0.500000000000000, 0.2703422189956631, 0.1024164803091, 0.01990973007891617, 0., 0., 0., 0., 0., 0., 0.};
		for(unsigned kx = 0; kx < nx; ++kx)
			REQUIRE(res[(nz/2)*ny*nz + (ny/2)*nx + kx] == Approx(row[kx]));
		// Check symmetry
		for(unsigned kz = 0; kz < nz; ++kz)
			for(unsigned ky = 0; ky < ny; ++ky)
				for(unsigned kx = 0; kx < nx/2; ++kx)
					REQUIRE(fabs(res[kz*ny*nz + ky*nx + kx] - 0.5) == Approx(fabs(res[kz*ny*nz + ky*nx + nx - kx - 1] - 0.5)));
		for(unsigned kz = 0; kz < nz; ++kz)
			for(unsigned ky = 0; ky < ny/2; ++ky)
				for(unsigned kx = 0; kx < nx; ++kx)
					REQUIRE(res[kz*ny*nz + ky*nx + kx] == Approx(res[kz*ny*nz + (ny - ky - 1)*nx + kx]));
		for(unsigned kz = 0; kz < nz/2; ++kz)
			for(unsigned ky = 0; ky < ny; ++ky)
				for(unsigned kx = 0; kx < nx; ++kx)
					REQUIRE(res[kz*ny*nz + ky*nx + kx] == Approx(res[(nz - kz - 1)*ny*nz + ky*nx + kx]));
		// Test single points
		std::vector<Point_3_base> ptVec = {Point_3_base(0., 0.5, 0.5), Point_3_base(0.5, 0.5, 0.5), 
			Point_3_base(1., 0.5, 0.5), Point_3_base(2., 2., 2.)};
		row = {1., 0.5, 0., 0.}; // Check values
		for(std::size_t k = 0; k < ptVec.size(); ++k)
			REQUIRE(testFilt(ptVec[k]) == Approx(row[k]));
		// Retest with the filter derivative (same results can be computed since it's a linear filter)
		SparseMatrix diffMat = testFilt.diffFilter(ptVec, rad);
		res = filterWithDiff(xVec, ptVec.size(), diffMat);
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(row[k]));
	}
}

TEST_CASE("Testing Filter3D with TOMesh tet input","[Filter3D]")
{
	// Read mesh from file
	std::ifstream meshFile("filter3dmesh.txt");
	REQUIRE(meshFile.good());
	// Read points
	std::size_t npts;
	meshFile >> npts;
	std::vector<Point_3_base> ptVec(npts);
	for(auto pit = ptVec.begin(); pit != ptVec.end(); ++pit)
	{
		double x, y, z;
		meshFile >> x >> y >> z;
		*pit = Point_3_base(x, y, z);
	}
	// Read element connectivity
	std::size_t nelems;
	meshFile >> nelems;
	std::vector<std::vector<std::size_t>> connectivity(nelems);
	for(auto cit = connectivity.begin(); cit != connectivity.end(); ++cit)
	{
		std::size_t n1, n2, n3, n4;
		meshFile >> n1 >> n2 >> n3 >> n4;
		*cit = {n1, n2, n3, n4};
	}
	// Read previously generated results
	std::vector<double> centroidRes(nelems);
	std::size_t ncentres;
	meshFile >> ncentres;
	REQUIRE(ncentres == nelems);
	for(auto cit = centroidRes.begin(); cit != centroidRes.end(); ++cit)
		meshFile >> *cit;
	std::vector<double> nodalRes(npts);
	std::size_t nnoderes;
	meshFile >> nnoderes;
	REQUIRE(nnoderes == npts);
	for(auto it = nodalRes.begin(); it != nodalRes.end(); ++it)
		meshFile >> *it;
	// Form mesh and test
	TOMesh3D testMesh(ptVec, connectivity);
	bool ptsAtCenter = true;
	SECTION("Testing centroid points")
	{
		double rad = 0.2;
		Filter3D<> testFilt0(&testMesh, rad);
		// Use assignment operator
		Filter3D<> testFilt = testFilt0;
		std::vector<double> xVec(testMesh.getNumElements());
		std::vector<Point_3_base> centroidPts(xVec.size());
		for(std::size_t k = 0; k < xVec.size(); ++k)
		{
			Point_3_base centroid = TOMeshProcessing::getElementCentroid3D(k, &testMesh);
			centroidPts[k] = centroid;
			if(centroid.y() < 0.5)
				xVec[k] = 0.;
			else
				xVec[k] = 1.;
		}
		std::vector<double> res = testFilt(xVec, rad);
		REQUIRE(res.size() == xVec.size());
		REQUIRE(res.size() == centroidRes.size());
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(centroidRes[k]));
		// Test outside region
		REQUIRE(testFilt(Point_3_base(2.,2.,2)) == Approx(0.));
		REQUIRE(testFilt(Point_3_base(-2.,2.,2.)) == Approx(0.));
		REQUIRE(testFilt(Point_3_base(2.,-2.,2.)) == Approx(0.));
		REQUIRE(testFilt(Point_3_base(2.,2.,-2.)) == Approx(0.));
		// Retest with the filter derivative (same results can be computed since it's a linear filter)
		SparseMatrix diffMat = testFilt.diffFilter(centroidPts, rad);
		res = filterWithDiff(xVec, centroidRes.size(), diffMat);
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(centroidRes[k]));
	}
	SECTION("Testing nodal points")
	{
		double rad = 0.2;
		Filter3D<> testFilt0(&testMesh, rad, false);
		// Use copy ctor
		Filter3D<> testFilt(testFilt0);
		std::vector<double> xVec(testMesh.getNumNodes());
		for(std::size_t k = 0; k < xVec.size(); ++k)
		{
			Point_3_base curNode = testMesh.getNode3D(k);
			if(curNode.y() < 0.5)
				xVec[k] = 0.;
			else
				xVec[k] = 1.;
		}
		std::vector<double> res = testFilt(xVec, rad);
		REQUIRE(res.size() == xVec.size());
		REQUIRE(res.size() == nodalRes.size());
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(nodalRes[k]));
		// Test outside region
		REQUIRE(testFilt(Point_3_base(2.,2.,2)) == Approx(0.));
		REQUIRE(testFilt(Point_3_base(-2.,2.,2.)) == Approx(0.));
		REQUIRE(testFilt(Point_3_base(2.,-2.,2.)) == Approx(0.));
		REQUIRE(testFilt(Point_3_base(2.,2.,-2.)) == Approx(0.));
		// Retest with the filter derivative (same results can be computed since it's a linear filter)
		SparseMatrix diffMat = testFilt.diffFilter(ptVec, rad);
		res = filterWithDiff(xVec, ptVec.size(), diffMat);
		for(std::size_t k = 0; k < res.size(); ++k)
			REQUIRE(res[k] == Approx(nodalRes[k]));
	}
}

