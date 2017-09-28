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

#include "postprocess.h"
#include "topoptrep.h"
#include "geometrytranslation.h"
#include "outputwriter.h"
#include "tomesh.h"
#include "helper.h"
#include "tomeshprocessing.h"
#include <vector>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Bbox_2.h>

namespace Topologies{
namespace PostProcess
{
	namespace
	{
		void replaceLargestPolygon(const std::vector<std::vector<Mesh_Segment_2> >& segVec1, 
															std::vector<std::vector<Mesh_Segment_2> >& segVec2)
		{
			// Find largest in 1
			double area = 0.;
			std::size_t klarge1 = 0;
			for(std::size_t k = 0; k < segVec1.size(); ++k)
			{
				double curArea = GeometryTranslation::computeSignedArea(segVec1[k]);
				if(curArea > area)
				{
					area = curArea;
					klarge1 = k;
				}
			}
			area = 0.;
			std::size_t klarge2 = 0;
			for(std::size_t k = 0; k < segVec2.size(); ++k)
      {
				double curArea = GeometryTranslation::computeSignedArea(segVec2[k]);
				if(curArea > area)
				{
					area = curArea;
					klarge2 = k;
				}
      }
			if(!segVec1.empty() && !segVec2.empty())
				segVec2[klarge2] = segVec1[klarge1];
		}
	
		bool isSegOn(const Mesh_Segment_2& inSeg, const std::vector<Mesh_Segment_2>& segVec, double tol2)
		{
			for(std::size_t k = 0; k < segVec.size(); ++k)
			{
				if(CGAL::squared_distance(segVec[k], inSeg.source()) < tol2 || 
					 CGAL::squared_distance(segVec[k], inSeg.target()) < tol2)
					return true;
			}
			return false;
		}

	} // Anonymous namespace for helper functions

	VoxelInterp::VoxelInterp(const std::vector<double>& inVA, unsigned inNx, unsigned inNy, unsigned inNz,
                double inW, double inL, double inH, double inThreshold) :
		voxelArray(inVA),
		nx(inNx), ny(inNy), nz(inNz),
		width(inW), length(inL), height(inH),
		threshold(inThreshold)
	{
		assert(voxelArray.size() == (nx+1)*(ny+1)*(nz+1));
	}

	Tr_GT::FT VoxelInterp::operator()(Tr_GT::Point_3 p) const
	{
		if(p.x() < 0. || p.y() < 0. || p.z() < 0.)
			return Tr_GT::FT(threshold);
		if(p.x() > width || p.y() > length || p.z() > height)
			return Tr_GT::FT(threshold);
		// Trilinear interpolation
		double dx = width/(double)nx, dy = length/(double)ny, dz = height/(double)nz;
		unsigned kpx = (unsigned)(p.x()/dx), kpy = (unsigned)(p.y()/dy), kpz = (unsigned)(p.z()/dz);
		if(kpx == nx)
			--kpx;
		if(kpy == ny)
			--kpy;
		if(kpz == nz)
			--kpz;
		double interpVal = 0.;
		for(unsigned kx = 0; kx < 2; ++kx)
		{
			double curxe = dx*(double)(kpx + kx);
			double xval = 0.;
			for(unsigned ky = 0; ky < 2; ++ky)
			{
				double curye = dy*(double)(kpy + ky);
				double yval = 0.;
				for(unsigned kz = 0; kz < 2; ++kz)
				{
					double curze = dz*(double)(kpz + kz);
					std::size_t kpxe = kpx + kx, kpye = kpy + ky, kpze = kpz + kz;
					std::size_t kp = kpze*(nx + 1)*(ny + 1) + kpye*(nx + 1) + kpxe;
					double zval = voxelArray[kp];
					double zdist = 1. - fabs(curze - p.z())/dz;
					yval += zdist*zval;
				}
				double ydist = 1. - fabs(curye - p.y())/dy;
				xval += ydist*yval;
			}
			double xdist = 1. - fabs(curxe - p.x())/dx;
			interpVal += xdist*xval;
		}
		return Tr_GT::FT(threshold - interpVal);
/*
		if(interpVal > threshold)
			return Tr_GT::FT(-1.);
		else if(interpVal < threshold)
			return Tr_GT::FT(1.);
		return Tr_GT::FT(0.); // interpVal == threshold*/
	}

	namespace
	{
		std::vector<Point_3_base> getTestPoints(const TOMesh* const inMesh, std::size_t ke)
		{
			std::vector<std::size_t> curElem = inMesh->getElementConnectivity(ke);
			if(curElem.size() == 4) // Tet
				return {TOMeshProcessing::getElementCentroid3D(ke, inMesh)};
			// Otherwise, compute Delaunay tesselation of element and use centroids of those tets
			std::vector<Point_3_base> ptVec(curElem.size());
			for(std::size_t k = 0; k < ptVec.size(); ++k)
				ptVec[k] = inMesh->getNode3D(curElem[k]);
			CGAL::Triangulation_3<Mesh_K> dtess(ptVec.begin(), ptVec.end());
			ptVec.clear();
			ptVec.reserve(5*dtess.number_of_finite_cells());
			for(CGAL::Triangulation_3<Mesh_K>::Finite_cells_iterator cit = dtess.finite_cells_begin();
          cit != dtess.finite_cells_end(); ++cit)
      {
				std::vector<Point_3_base> verts = {cit->vertex(0)->point(), cit->vertex(1)->point(), 
					cit->vertex(2)->point(), cit->vertex(3)->point()};
				// Centroid
				ptVec.push_back(CGAL::centroid(verts[0], verts[1], verts[2], verts[3]));
				// Points close to vertices
				for(unsigned k = 0; k < 4; ++k)
				{
					Mesh_K::Vector_3 rxi = verts[(k+1)%4] - verts[k], reta = verts[(k+2)%4] - verts[k], rzeta = verts[(k+3)%4] - verts[k];
					double ac = 0.05;
					ptVec.push_back(verts[k] + ac*(rxi + reta + rzeta)); // Point near vertex
					ac = 0.3;
					ptVec.push_back(verts[k] + ac*(rxi + reta + rzeta)); // Point near face
				}
			}
			return ptVec;
		}
	}

	TetMeshInterp::TetMeshInterp(const std::vector<double>& inVA, const TOMesh* const inMesh, double inThreshold) :
		nodalValArray(inVA), 
		threshold(inThreshold),
		valid(false)
	{
		assert(nodalValArray.size() == inMesh->getNumNodes());
		std::vector<std::pair<TMI_Delaunay::Point, std::size_t>> ptVec(inMesh->getNumNodes());
		for(std::size_t k = 0; k < ptVec.size(); ++k)
			ptVec[k] = std::make_pair(inMesh->getNode3D(k), k);
		tess = TMI_Delaunay(ptVec.begin(), ptVec.end());
		if(tess.is_valid() && tess.number_of_vertices() >= 4 && tess.number_of_cells() >= 1)
		{
			valid = true;
			// Store delaunay cells that are in the original mesh
			for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			{
				std::vector<Point_3_base> tmpTestVec = getTestPoints(inMesh, k);
				testPtVec.insert(testPtVec.end(), tmpTestVec.begin(), tmpTestVec.end());
			}
			generateCellSet();
		}
		if(tess.number_of_finite_cells() != cellSet.size())
		{
			std::cout << "Warning: Didn't find all Delaunay tets in original mesh, hopefully it's concave!" << std::endl;
			std::cout << "Num delaunay: " << tess.number_of_finite_cells() << ", num cell set: " << cellSet.size() << std::endl;
		}
	}

	TetMeshInterp::TetMeshInterp(const TetMeshInterp& copy) :
		threshold(copy.threshold),
		nodalValArray(copy.nodalValArray),
		tess(copy.tess),
		testPtVec(copy.testPtVec),
		valid(copy.valid)
	{
		if(valid)
			generateCellSet();
	}

	TetMeshInterp::TetMeshInterp(TetMeshInterp&& copy)
	{
		swap(copy);
	}

	TetMeshInterp& TetMeshInterp::operator=(TetMeshInterp rhs)
	{
		swap(rhs);
		return *this;
	}

	void TetMeshInterp::swap(TetMeshInterp& arg2)
	{
		std::swap(threshold, arg2.threshold);
		nodalValArray.swap(arg2.nodalValArray);
		tess.swap(arg2.tess);
		testPtVec.swap(arg2.testPtVec);
		cellSet.swap(arg2.cellSet);
		std::swap(valid, arg2.valid);
	}
	
	void TetMeshInterp::generateCellSet()
	{
		cellSet.clear();
		for(auto it = testPtVec.begin(); it != testPtVec.end(); ++it)
		{
			TMI_Delaunay::Cell_handle ch = tess.locate(*it);
			cellSet.insert(ch);
		}
	}
	Tr_GT::FT TetMeshInterp::operator()(Tr_GT::Point_3 p) const
	{
		if(!valid)
			return 1.;
		// Find enclosing cell
		TMI_Delaunay::Point pbP(p.x(), p.y(), p.z());
		TMI_Delaunay::Cell_handle ch = tess.locate(pbP);
		// Check whether or not it was contained in the original mesh
		if(tess.is_infinite(ch) || (cellSet.find(ch) == cellSet.end()))
			return 1.; // Outside region
		// Compute area coordinates and get interpolated value
		Point_3_base ac = computeAreaCoords(ch, p);
		unsigned vid = ch->vertex(0)->info();
		double interpVal = (1. - ac.x() - ac.y() - ac.z())*nodalValArray[vid];
		for(unsigned k = 1; k < 4; ++k)
		{
			vid = ch->vertex(k)->info();
			interpVal += ac[k - 1]*nodalValArray[vid];
		}
		return Tr_GT::FT(threshold - interpVal);
/*
		if(interpVal > threshold)
			return Tr_GT::FT(-1.);
		else if(interpVal < threshold)
			return Tr_GT::FT(1.);
		return Tr_GT::FT(0.); // interpVal == threshold*/
	}

	Point_3_base TetMeshInterp::computeAreaCoords(const TMI_Delaunay::Cell_handle ch, const Point_3_base& p) const
	{
		Point_3_base p000 = ch->vertex(0)->point(), p100 = ch->vertex(1)->point(), 
					 p010 = ch->vertex(2)->point(), p001 = ch->vertex(3)->point();
		Mesh_K::Vector_3 rxi = p100 - p000, reta = p010 - p000, rzeta = p001 - p000;
		CGAL::Aff_transformation_3<Mesh_K> atran(rxi.x(), reta.x(), rzeta.x(), p000.x(), 
			rxi.y(), reta.y(), rzeta.y(), p000.y(), 
			rxi.z(), reta.z(), rzeta.z(), p000.z());
		CGAL::Aff_transformation_3<Mesh_K> invTrans(atran.inverse());
		return invTrans(p);
	}

	std::vector<Mesh_Segment_2> meshSegsIsoSurf2d(const std::vector<double>& valVec, const std::vector<Point_2_base>& ptVec, double threshold)
	{
		// Compute tolerance based on max distance between points
		// This sets the minimum segment length
		CGAL::Bbox_2 bb = CGAL::bbox_2(ptVec.begin(), ptVec.end());
		double xd = bb.xmax() - bb.xmin(), yd = bb.ymax() - bb.ymin();
		double tol = sqrt(xd*xd + yd*yd)*1e-10;
		// Form Delaunay triangulation of 3d points
		typedef CGAL::Projection_traits_xy_3<Mesh_K>  Ptxy3;
		typedef CGAL::Delaunay_triangulation_2<Ptxy3> Delaunay;
		std::vector<Point_3_base> terrainPtVec(ptVec.size());
		if(valVec.size() != ptVec.size())
		{
			std::cout << "Error in meshSegsIsoSurf2d: valVec size (" << valVec.size() << ") doesn't match ptVec size (" << ptVec.size() << ")" << std::endl;
			std::vector<Mesh_Segment_2> bsVec;
			return bsVec;
		}
		for(std::size_t k = 0; k < terrainPtVec.size(); ++k)
			terrainPtVec[k] = Point_3_base(ptVec[k].x(), ptVec[k].y(), valVec[k]);
		Delaunay dt(terrainPtVec.begin(), terrainPtVec.end());
		// Now find intersections of triangles with threshold plane
		CGAL::Plane_3<Mesh_K>::Plane_3 threshPlane(0., 0., 1., -threshold);
		std::vector<Mesh_Segment_2> outVec;
		for(Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit)
		{
			CGAL::Triangle_3<Mesh_K> curTri = dt.triangle(fit);
			auto res = CGAL::intersection(threshPlane, curTri);
			if(res)
			{
				const Mesh_K::Segment_3* sb = boost::get<Mesh_K::Segment_3>(&*res);
				if(sb)
				{
					// Found intersection!
					Point_3_base s(sb->source()), t(sb->target());
					// Check length
					Mesh_Segment_2 newSeg(Point_2_base(s.x(), s.y()), Point_2_base(t.x(), t.y()));
					if(sqrt(newSeg.squared_length()) > tol)
						outVec.push_back(Mesh_Segment_2(Point_2_base(s.x(), s.y()), Point_2_base(t.x(), t.y())));
				}
			}
		}
		return outVec;
	}

	void linearExtrudeIsoSurf(const TopOptRep& inTORep, const std::string& fileName, double height)
	{
		std::vector<std::size_t> sizes;
		std::vector<double> realVec;
		inTORep.getDataSize(sizes);
		inTORep.getRealRep(realVec);
		if(sizes.size() != 2)
			return;
		typedef CGAL::Implicit_surface_3<Tr_GT, VoxelInterp> Voxel_Surface_3;
		SM3_Tr tr;
		C2t3 c2t3(tr);
		std::size_t nx = sizes[0], ny = sizes[1], nz = 1;
		double width = 1., length = 1.;
		VoxelInterp voxi(realVec, nx, ny, nz, width, length, height, 0.5);
		Tr_GT::Point_3 bounding_sphere_center(width*0.5, length*0.5, height*0.5);
		Tr_GT::FT bounding_sphere_squared_radius = 1.2*(width*width*0.25 + length*length*0.25 + height*height*0.25);
		Tr_GT::Sphere_3 bounding_sphere(bounding_sphere_center, bounding_sphere_squared_radius);
		Voxel_Surface_3 surface(voxi, bounding_sphere, 1e-5);
		double dx = width/(double)nx;
		CGAL::Surface_mesh_default_criteria_3<SM3_Tr> criteria(30., 0.5*dx, 0.5*dx);
		CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
		OutputWriter::plotSTL(c2t3, fileName);
	}

	void linearExtrude(const TopOptRep& inTORep, const std::string& fileName, double length)
	{
		// Get 2d segments and order them
		std::vector<Mesh_Segment_2> segVec;
		inTORep.get2DSegments(segVec);
		std::vector<std::vector<Mesh_Segment_2>> orderedSegVec;
		GeometryTranslation::orderMeshSegments(segVec, orderedSegVec);
		std::cout << "linearExtrude: Num components: " << orderedSegVec.size() << std::endl;
		// Simplify and subdivide to smooth
		std::vector<std::vector<Mesh_Segment_2>> subdSegVec(orderedSegVec.size());
		std::vector<Mesh_Segment_2> boundaryVec;
		inTORep.getBoundary(boundaryVec);
		for(std::size_t k = 0; k < orderedSegVec.size(); ++k)
		{
			simplify2d(orderedSegVec[k], subdSegVec[k]);
			subdivision2d(subdSegVec[k], orderedSegVec[k], boundaryVec);
		}
		// Mesh new result
		std::vector<GenericMaterial> matVec(1); // empty
		GeometryTranslation::MesherData meshParams;
		meshParams.triMeshEdgeAngle = 10.;
		meshParams.triMeshEdgeSize = 0.1; // Edge size doesn't really matter
		CDT_2 newMesh = GeometryTranslation::mesh2D(orderedSegVec, matVec, meshParams);
		// Finally extrude
		OutputWriter::plotSTL(GeometryTranslation::linearExtrude(newMesh, length), fileName, true);
	}

	void simplify2d(const std::vector<Mesh_Segment_2>& inVec, std::vector<Mesh_Segment_2>& outVec)
	{
		outVec.clear();
		if(inVec.empty())
			return;
		bool isCycle = TOL_EQ(inVec.front().source().x(), inVec.back().target().x(), 1e-12) &&
				TOL_EQ(inVec.front().source().y(), inVec.back().target().y(), 1e-12);
		Point_2_base curSource = inVec[0].source();
		bool lastCheck;
		std::size_t decrement = isCycle ? 0 : 1; 
		for(std::size_t k = 0; k < inVec.size() - decrement; ++k)
		{
			std::size_t kp1 = (k + 1) % inVec.size();
			double dotp = inVec[k].to_vector() * inVec[kp1].to_vector();
			double l1 = sqrt(inVec[k].squared_length()), l2 = sqrt(inVec[kp1].squared_length());
			lastCheck = !TOL_EQ(dotp/(l1*l2), 1., 1e-12);
			if(lastCheck)
			{
				// Construct new segment
				Mesh_Segment_2 newSeg(curSource, inVec[k].target());
				outVec.push_back(newSeg);
				curSource = inVec[k].target();
			}
		}
		if(isCycle && (!lastCheck && !outVec.empty())) // Replace first segment's source point
			outVec[0] = Mesh_Segment_2(curSource, outVec[0].target());
		else if(!isCycle) // inVec is 1 straight line
			outVec.push_back(Mesh_Segment_2(curSource, inVec.back().target()));
	}

	void subdivision2d(const std::vector<Mesh_Segment_2>& coarse, std::vector<Mesh_Segment_2>& fine, 
										const std::vector<Mesh_Segment_2>& boundaryVec)
	{
		fine.clear();
		double t = 0.25;
		std::size_t ncp = coarse.size();
		std::size_t nignored = 0;
		double tol = 1e-10;
		for(std::size_t k = 0; k < ncp; ++k)
		{
			if(!isSegOn(coarse[k], boundaryVec, tol*tol))
			{
				std::size_t km1 = (k + (ncp - 1)) % ncp;
				Point_2_base ep0;
				if(!isSegOn(coarse[km1], boundaryVec, tol*tol))
					ep0 = coarse[k].source() + t*coarse[k].to_vector();
				else
					ep0 = coarse[k].source();
				Point_2_base ep1 = coarse[k].source() + (1. - t)*coarse[k].to_vector();
				fine.push_back(Mesh_Segment_2(ep0, ep1));
				std::size_t kp1 = (k + 1) % ncp;
				ep0 = ep1;
				if(!isSegOn(coarse[kp1], boundaryVec, tol*tol))
					ep1 = coarse[kp1].source() + t*coarse[kp1].to_vector();
				else
					ep1 = coarse[kp1].source();
				fine.push_back(Mesh_Segment_2(ep0, ep1));
			}
			else
			{
				++nignored;
				fine.push_back(coarse[k]);
			}
		}
	}
}
}

