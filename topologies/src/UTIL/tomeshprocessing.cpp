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

#include <vector>
#include "tomesh.h"
#include <CGAL/Triangulation_3.h>

namespace Topologies{
namespace TOMeshProcessing
{
	namespace
	{
		double computeTriOrQuadArea(const std::vector<std::size_t>& elemConn, const TOMesh* const inMesh)
		{
			// Uses CGAL::Polygon_2 to compute area and ensure correct node ordering of quads
			std::vector<Mesh_K::Point_2> tmpPtVec(elemConn.size());
			for(std::size_t kn = 0; kn < elemConn.size(); ++kn)
				tmpPtVec[kn] = inMesh->getNode2D(elemConn[kn]);
			CGAL::Polygon_2<Mesh_K> tmpPoly(tmpPtVec.begin(), tmpPtVec.end());
			std::size_t kn = 1;
			while(!tmpPoly.is_simple() && kn < tmpPtVec.size())
			{
				// Out of order nodes, assumed this is a quadrilateral so one swap will order the nodes correctly
				std::swap(tmpPtVec[0], tmpPtVec[kn]);
				tmpPoly = CGAL::Polygon_2<Mesh_K>(tmpPtVec.begin(), tmpPtVec.end());
				// Swap back
				std::swap(tmpPtVec[0], tmpPtVec[kn++]);
			}
			return fabs(tmpPoly.area());
		}

		double computeTetVolume(const std::vector<std::size_t>& elemConn, const TOMesh* const inMesh)
		{
			return fabs(CGAL::volume(inMesh->getNode3D(elemConn[0]), inMesh->getNode3D(elemConn[1]), 
															inMesh->getNode3D(elemConn[2]), inMesh->getNode3D(elemConn[3])));
		}

		double computeHexVolume(const std::vector<std::size_t>& elemConn, const TOMesh* const inMesh)
		{
			// Requires some pre-processing to ensure correct node ordering
			// A 3d Delaunay tesselation is computed, then the tetrahedra volumes computed and summed
			std::vector<Point_3_base> ptVec(elemConn.size());
			for(std::size_t k = 0; k < ptVec.size(); ++k)
				ptVec[k] = inMesh->getNode3D(elemConn[k]);
			CGAL::Triangulation_3<Mesh_K> dtess(ptVec.begin(), ptVec.end());
			assert(dtess.is_valid()); // Check validity of tesselation
			// Sum tetrahedra volumes
			double volsum = 0.;
			for(CGAL::Triangulation_3<Mesh_K>::Finite_cells_iterator cit = dtess.finite_cells_begin();
					cit != dtess.finite_cells_end(); ++cit)
			{
				volsum += fabs(CGAL::volume(cit->vertex(0)->point(), cit->vertex(1)->point(), 
																		cit->vertex(2)->point(), cit->vertex(3)->point()));
			}
			return volsum;
		}
	}

  double computeElementVolume(std::size_t elemID, const TOMesh* const inMesh)
	{
		std::vector<std::size_t> curElem = inMesh->getElementConnectivity(elemID);
		if(inMesh->dimNum() == 2 && (curElem.size() == 3 || curElem.size() == 4))
			return computeTriOrQuadArea(curElem, inMesh);
		else if(inMesh->dimNum() == 3)
		{
			if(curElem.size() == 4)
				return computeTetVolume(curElem, inMesh);
			else if(curElem.size() == 8)
				return computeHexVolume(curElem, inMesh);
		}
		return 0.;
	}

	Point_2_base getElementCentroid2D(std::size_t elemID, const TOMesh* const inMesh)
  {
		std::vector<std::size_t> curElem = inMesh->getElementConnectivity(elemID);
		assert(curElem.size() == 3 || curElem.size() == 4); // Must be triangle or quad
		Point_2_base p1 = inMesh->getNode2D(curElem[0]);
		Point_2_base p2 = inMesh->getNode2D(curElem[1]);
		Point_2_base p3 = inMesh->getNode2D(curElem[2]);
		if(curElem.size() == 3) // Triangular element
			return CGAL::centroid(p1,p2,p3);
		// Quadrilateral element
		Point_2_base p4 = inMesh->getNode2D(curElem[3]);
		return CGAL::centroid(p1,p2,p3,p4);
  }

	Point_3_base getElementCentroid3D(std::size_t elemID, const TOMesh* const inMesh)
	{
		std::vector<std::size_t> connectVec = inMesh->getElementConnectivity(elemID);
		Mesh_K::Vector_3 sumv(0., 0., 0.);
		for(std::size_t kp = 0; kp < connectVec.size(); ++kp)
			sumv = sumv + (inMesh->getNode3D(connectVec[kp]) - CGAL::ORIGIN);
		return CGAL::ORIGIN + sumv/(double)connectVec.size();
	}
}
}
