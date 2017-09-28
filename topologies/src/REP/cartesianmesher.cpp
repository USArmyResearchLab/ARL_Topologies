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

#include "cartesianmesher.h"

namespace Topologies{
namespace CartesianMesher
{
	namespace
	{
		std::vector<Point_2_base> getPointVec(unsigned nx, unsigned ny, double width, double height)
		{
			std::vector<Point_2_base> ptVec((nx + 1)*(ny + 1));
			double dx = width/(double)nx, dy = height/(double)ny;
			for(std::size_t ky = 0; ky <= ny; ++ky)
			{
				double cury = (double)ky*dy;
				for(std::size_t kx = 0; kx <= nx; ++kx)
				{
					double curx = (double)kx*dx;
					ptVec[ky*(nx + 1) + kx] = Point_2_base(curx, cury);
				}
			}
			return ptVec;
		}

		std::vector<Point_3_base> getPointVec(unsigned nx, unsigned ny, unsigned nz, double width, double length, double height)
		{
			std::vector<Point_3_base> ptVec((nx+1)*(ny+1)*(nz+1));
			double dx = width/(double)nx, dy = length/(double)ny, dz = height/(double)nz;
			for(std::size_t kz = 0; kz <= nz; ++kz)
			{
				double curz = (double)kz*dz;
				for(std::size_t ky = 0; ky <= ny; ++ky)
				{
					double cury = (double)ky*dy;
					for(std::size_t kx = 0; kx <= nx; ++kx)
					{
						double curx = (double)kx*dx;
						ptVec[kz*(ny+1)*(nx+1) + ky*(nx+1) + kx] = Point_3_base(curx, cury, curz);
					}
				}
			}
			return ptVec;
		}

		std::unique_ptr<TOMesh2D> getQuadMesh(unsigned nx, unsigned ny, double width, double height, const std::vector<double>& optVals)
		{
			// Set up element connectivity
			std::vector<std::vector<std::size_t>> elemConnVec(nx*ny);
			for(std::size_t ky = 0; ky < ny; ++ky)
			{
				for(std::size_t kx = 0; kx < nx; ++kx)
				{
					std::vector<std::size_t>& curElem = elemConnVec[ky*nx + kx];
					curElem.resize(4);
					curElem[0] = ky*(nx + 1) + kx; // Bottom left corner
					curElem[1] = ky*(nx + 1) + kx + 1; // Bottom right
					curElem[2] = (ky + 1)*(nx + 1) + kx + 1; // Top right
					curElem[3] = (ky + 1)*(nx + 1) + kx; // Top left
				}
			}
			std::vector<Point_2_base> ptVec = getPointVec(nx, ny, width, height);
			return std::unique_ptr<TOMesh2D>(new TOMesh2D(ptVec, elemConnVec, optVals));
		}

		std::unique_ptr<TOMesh2D> getTriMesh(unsigned nx, unsigned ny, double width, double height, const std::vector<double>& optVals)
		{
			// Set up element connectivity
			std::vector<std::vector<std::size_t>> elemConnVec(2*nx*ny);
			std::size_t k = 0;
			for(std::size_t ky = 0; ky < ny; ++ky)
			{
				for(std::size_t kx = 0; kx < nx; ++kx)
				{
					std::vector<std::size_t>& curElem1 = elemConnVec[k++];
					curElem1.resize(3);
					curElem1[0] = ky*(nx + 1) + kx;
					curElem1[1] = ky*(nx + 1) + kx + 1;
					curElem1[2] = (ky + 1)*(nx + 1) + kx + 1;
					std::vector<std::size_t>& curElem2 = elemConnVec[k++];
					curElem2.resize(3);
					curElem2[0] = ky*(nx + 1) + kx;
					curElem2[1] = (ky + 1)*(nx + 1) + kx + 1;
					curElem2[2] = (ky + 1)*(nx + 1) + kx;
				}
			}
			std::vector<Point_2_base> ptVec = getPointVec(nx, ny, width, height);
			// Double optVal vec
			std::vector<double> elemValVec(optVals.size()*2);
			for(std::size_t k = 0; k < optVals.size(); ++k)
			{
				elemValVec[2*k] = optVals[k];
				elemValVec[2*k + 1] = optVals[k];
			}
			return std::unique_ptr<TOMesh2D>(new TOMesh2D(ptVec, elemConnVec, elemValVec));
		}

		std::unique_ptr<TOMesh3D> getHexMesh(unsigned nx, unsigned ny, unsigned nz, double width, double length, double height, const std::vector<double>& optVals)
		{
			// Set up element connectivity
			std::vector<std::vector<std::size_t>> elemConnVec(nx*ny*nz);
			for(std::size_t kz = 0; kz < nz; ++kz)
			{
				for(std::size_t ky = 0; ky < ny; ++ky)
				{
					for(std::size_t kx = 0; kx < nx; ++kx)
					{
						std::vector<std::size_t>& curElem = elemConnVec[kz*ny*nx + ky*nx + kx];
						curElem.resize(8);
						curElem[0] = kz*(ny+1)*(nx+1) + ky*(nx+1) + kx;
						curElem[1] = kz*(ny+1)*(nx+1) + ky*(nx+1) + kx + 1;
						curElem[2] = kz*(ny+1)*(nx+1) + (ky+1)*(nx+1) + kx + 1;
						curElem[3] = kz*(ny+1)*(nx+1) + (ky+1)*(nx+1) + kx;
						curElem[4] = (kz+1)*(ny+1)*(nx+1) + ky*(nx+1) + kx;
						curElem[5] = (kz+1)*(ny+1)*(nx+1) + ky*(nx+1) + kx + 1;
						curElem[6] = (kz+1)*(ny+1)*(nx+1) + (ky+1)*(nx+1) + kx + 1;
						curElem[7] = (kz+1)*(ny+1)*(nx+1) + (ky+1)*(nx+1) + kx;
					}
				}
			}
			std::vector<Point_3_base> ptVec = getPointVec(nx, ny, nz, width, length, height);
			return std::unique_ptr<TOMesh3D>(new TOMesh3D(ptVec, elemConnVec, optVals));
		}

		std::unique_ptr<TOMesh3D> getTetMesh(unsigned nx, unsigned ny, unsigned nz, double width, double length, double height, const std::vector<double>& optVals)
		{
			// Set up element connectivity
			std::vector<std::vector<std::size_t>> elemConnVec(6*nx*ny*nz);
			std::size_t k = 0;
			for(std::size_t kz = 0; kz < nz; ++kz)
			{
				for(std::size_t ky = 0; ky < ny; ++ky)
				{
					for(std::size_t kx = 0; kx < nx; ++kx)
					{
						std::size_t k0 = kz*(ny+1)*(nx+1) + ky*(nx+1) + kx,
										 k1 = kz*(ny+1)*(nx+1) + ky*(nx+1) + kx + 1,
										 k2 = kz*(ny+1)*(nx+1) + (ky+1)*(nx+1) + kx,
										 k3 = kz*(ny+1)*(nx+1) + (ky+1)*(nx+1) + kx + 1,
										 k4 = (kz+1)*(ny+1)*(nx+1) + ky*(nx+1) + kx,
										 k5 = (kz+1)*(ny+1)*(nx+1) + ky*(nx+1) + kx + 1,
										 k6 = (kz+1)*(ny+1)*(nx+1) + (ky+1)*(nx+1) + kx,
										 k7 = (kz+1)*(ny+1)*(nx+1) + (ky+1)*(nx+1) + kx + 1;
						// Index values come from the Delaunay tesselation of a unit cube
						std::vector<std::size_t>& curElem1 = elemConnVec[k++];
						curElem1.resize(4);
						curElem1[0] = k4;
						curElem1[1] = k0;
						curElem1[2] = k2;
						curElem1[3] = k1;
						std::vector<std::size_t>& curElem2 = elemConnVec[k++];
						curElem2.resize(4);
						curElem2[0] = k6;
						curElem2[1] = k4;
						curElem2[2] = k2;
						curElem2[3] = k1;
						std::vector<std::size_t>& curElem3 = elemConnVec[k++];
						curElem3.resize(4);
						curElem3[0] = k6;
						curElem3[1] = k5;
						curElem3[2] = k4;
						curElem3[3] = k1;
						std::vector<std::size_t>& curElem4 = elemConnVec[k++];
						curElem4.resize(4);
						curElem4[0] = k6;
						curElem4[1] = k3;
						curElem4[2] = k5;
						curElem4[3] = k1;
						std::vector<std::size_t>& curElem5 = elemConnVec[k++];
						curElem5.resize(4);
						curElem5[0] = k6;
						curElem5[1] = k2;
						curElem5[2] = k3;
						curElem5[3] = k1;
						std::vector<std::size_t>& curElem6 = elemConnVec[k++];
						curElem6.resize(4);
						curElem6[0] = k6;
						curElem6[1] = k7;
						curElem6[2] = k5;
						curElem6[3] = k3;
					}
				}
			}
			// Sextuple optVal vec
			std::vector<double> elemValVec(optVals.size()*6);
			for(std::size_t k = 0; k < optVals.size(); ++k)
			{
				elemValVec[6*k + 0] = optVals[k];
				elemValVec[6*k + 1] = optVals[k];
				elemValVec[6*k + 2] = optVals[k];
				elemValVec[6*k + 3] = optVals[k];
				elemValVec[6*k + 4] = optVals[k];
				elemValVec[6*k + 5] = optVals[k];
			}
			std::vector<Point_3_base> ptVec = getPointVec(nx, ny, nz, width, length, height);
			return std::unique_ptr<TOMesh3D>(new TOMesh3D(ptVec, elemConnVec, elemValVec));
		}
	}

	std::unique_ptr<TOMesh2D> generateMesh(MeshElementType inMET, unsigned nx, unsigned ny, double width, double height,
			const std::vector<double>& optVals)
	{
		assert(optVals.size() == nx*ny);
		assert(inMET == metTri || inMET == metQuad);
		std::unique_ptr<TOMesh2D> outMesh;
		if(inMET == metTri)
			outMesh = getTriMesh(nx, ny, width, height, optVals);
		else if(inMET == metQuad)
			outMesh = getQuadMesh(nx, ny, width, height, optVals);
		return outMesh;
	}

	std::unique_ptr<TOMesh3D> generateMesh(MeshElementType inMET, unsigned nx, unsigned ny, unsigned nz,
			double width, double length, double height,
			const std::vector<double>& optVals)
	{
		assert(optVals.size() == nx*ny*nz);
		assert(inMET == metTet || inMET == metHex);
		std::unique_ptr<TOMesh3D> outMesh;
		if(inMET == metTet)
			outMesh = getTetMesh(nx, ny, nz, width, length, height, optVals);
		else if(inMET == metHex)
			outMesh = getHexMesh(nx, ny, nz, width, length, height, optVals);
		return outMesh;
	}
}
}
