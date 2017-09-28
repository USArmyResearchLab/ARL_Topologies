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

#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <string>
#include <vector>
#include <unordered_set>
#include <set>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include "cgal_types.h"

namespace Topologies{
class TopOptRep;
class TOMesh;

//! Classes and functions to aid in post-processing, such as interpolation, extrusion, subdivision, and simplification
/*! This namespace contains classes and functions to aid in post-processing.  Two classes are used as helpers for
 *  interfacing with CGAL.  VoxelInterp provides an interpolation scheme within a voxel for generating surfaces
 *  from VoxelRep.  Similarly, TetMeshInterp provides an interpolation scheme within a tetrahedron for 
 *  post-processing of TetMesh3D.  
 *
 *  In addition, there are routines for linear extrusion of 2D designs to 3D designs, 2d simplification
 *  (i.e. merging of colinear line segments), and 2d subdivision smoothing.  
 */
namespace PostProcess
{
	//! An interface class for generating surfaces using CGAL from VoxelRep topologies
	/*! The interpolation scheme is trilinear (https://en.wikipedia.org/wiki/Trilinear_interpolation) */
	class VoxelInterp
	{
		public:
		//! Constructor taking a vector of scalars denoting voxel values, the 3d size of the array, the 3d dimensions, and the threshold
		VoxelInterp(const std::vector<double>& inVA, unsigned inNx, unsigned inNy, unsigned inNz, 
								double inW, double inL, double inH, double inThreshold);
		//! Evaluation operator: Returns the interpolated value at an arbitrary point within the voxelized geometry
		Tr_GT::FT operator()(Tr_GT::Point_3 p) const;

		private:
		std::vector<double> voxelArray;
		double width, length, height; // Region physical size
		unsigned nx, ny, nz;
		double threshold;
	};

	//! An interface class for generating surfaces using CGAL from TetMesh3D topologies
	/*! The interpolation is linear, though remeshed using the tetrahedra centroids.  In other words,
		*  a tetrahedral mesh is passed in with scalar values attached to each element.  The centroids
		*  of the elements are computed and a new mesh is generated using a Delaunay tesselation.  The
		*  new mesh can then use the scalar value as nodal values, facilitating a linear interpolation scheme.
		*/
	class TetMeshInterp
	{
	public:
		//! Constructor taking a vector of scalars denoding tetrahedra values, a tetrahedra mesh, and a threshold
		TetMeshInterp(const std::vector<double>& inVA, const TOMesh* const inMesh, double inThreshold);
		TetMeshInterp(const TetMeshInterp& copy);
		TetMeshInterp(TetMeshInterp&& copy);
		TetMeshInterp& operator=(TetMeshInterp rhs);
		void swap(TetMeshInterp& arg2);
		//! Evaluation operator: Returns the interpolated value at an arbitrary point within the voxelized geometry
		Tr_GT::FT operator()(Tr_GT::Point_3 p) const;
	private:
		typedef CGAL::Triangulation_vertex_base_with_info_3<std::size_t, Mesh_K> TMI_Vb3;
		typedef CGAL::Triangulation_data_structure_3<TMI_Vb3> TMI_Tds;
		typedef CGAL::Delaunay_triangulation_3<Mesh_K, TMI_Tds, CGAL::Fast_location> TMI_Delaunay;
		Point_3_base computeAreaCoords(const TMI_Delaunay::Cell_handle ch, const Point_3_base& p) const;
		void generateCellSet();
		struct Cell_handle_hash
		{
			std::size_t operator()(TMI_Delaunay::Cell_handle const& s) const
			{
				std::size_t const h1 ( std::hash<unsigned>()(s->vertex(0)->info()) );
				std::size_t const h2 ( std::hash<unsigned>()(s->vertex(1)->info()) );
				std::size_t const h3 ( std::hash<unsigned>()(s->vertex(2)->info()) );
				std::size_t const h4 ( std::hash<unsigned>()(s->vertex(3)->info()) );
				return h1 ^ (h2 ^ (h3 ^ (h4 << 1) << 1) << 1);
			}
		};
		struct Cell_handle_equal
		{
			bool operator()(TMI_Delaunay::Cell_handle const& x, TMI_Delaunay::Cell_handle const& y) const
			{
				return (x->vertex(0)->info() == y->vertex(0)->info()) && 
					(x->vertex(1)->info() == y->vertex(1)->info()) &&
					(x->vertex(2)->info() == y->vertex(2)->info()) &&
					(x->vertex(3)->info() == y->vertex(3)->info());
			}
		};
		std::vector<double> nodalValArray;
		std::unordered_set<TMI_Delaunay::Cell_handle, Cell_handle_hash, Cell_handle_equal> cellSet;
		TMI_Delaunay tess;
		std::vector<Point_3_base> testPtVec;
		double threshold;
		bool valid;
	};
	//! A function to merge collinear line segments in a polygon
	/*! Input `inVec` is assumed to be a well-formed polygon, i.e. all segments are connected in order.
	 *  The output `outVec` consists of a simplified polygon, where collinear segments are merged as:
   *  .-----.------.----. becomes .-----------------. 
   */
	void simplify2d(const std::vector<Mesh_Segment_2>& inVec, std::vector<Mesh_Segment_2>& outVec);
	//! A function to perform 2D subdivision for smoothing polygons
	/*! A vector of boundary segments is also an input argument.  The function will not subdivide
   *  segments that are collinear with the boundary.
	 */
	void subdivision2d(const std::vector<Mesh_Segment_2>& coarse, std::vector<Mesh_Segment_2>& fine, 
										const std::vector<Mesh_Segment_2>& boundaryVec);
	//! A function to perform linear extrusion of a 2D design into a 3D geometry
	/*! This version uses the TopOptRep function get2DSegments and performs a straight-forward
	 *  linear extrusion of those segments.
	 */
	void linearExtrude(const TopOptRep& inTORep, const std::string& fileName, double length);
	//! A function to perform linear extrusion of a 2D design into a 3D geometry
	/*! This version constructs a voxel topology from the TopOptRep and uses CGAL's surface mesh
	 *  generation and the VoxelInterp methods.  The result will be smoother, but not an exact extrusion.
	 */
	void linearExtrudeIsoSurf(const TopOptRep& inTORep, const std::string& fileName, double height);
	//! Returns a set of line segments representing the 2d iso-surface of the height map defined by valVec and ptVec
	/*! This function constructs a 2D iso-surface using a height-map: the input arguments valVec and ptVec are
	 *  interpreted to be essentially a 3d point cloud, with valVec being the z-coordinate.  Input argument `threshold`
	 *  defines the z intercept of an x-y plane and the intersection of that plane and the height map is computed.
	 *  This intersection is then a set of 2D line segments, that become the iso-surface
	 */
	std::vector<Mesh_Segment_2> meshSegsIsoSurf2d(const std::vector<double>& valVec, const std::vector<Point_2_base>& ptVec, double threshold);
}
}
#endif

