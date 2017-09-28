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

#ifndef GEOMETRYTRANSLATION_H
#define GEOMETRYTRANSLATION_H

#include "cgal_types.h"
#include "topologiesdefs.h"
#include "genericmaterial.h"
#include "tomesh.h"
#include <vector>

namespace Topologies{
//! A collection of functions for generating meshes and processing geometry in general
/*! This namespace contains various functions for meshing simple regions such as
 *  pixel and voxel rects, as well as processing of 2d line segments. 
 */
namespace GeometryTranslation
{
	//! Simple struct to contain all mesher parameters
	struct MesherData
	{
		MesherData() : triMeshEdgeAngle(10.), triMeshEdgeSize(1.), 
			tetMeshEdgeSize(1.), tetMeshFacetAngle(10.), tetMeshFacetSize(1.), tetMeshFacetDistance(1.), 
			tetMeshCellRadiusEdgeRatio(1.), tetMeshCellSize(1.) {}
		Real triMeshEdgeAngle, triMeshEdgeSize; //2D
		Real tetMeshEdgeSize, tetMeshFacetAngle, tetMeshFacetSize, tetMeshFacetDistance,
			tetMeshCellRadiusEdgeRatio, tetMeshCellSize; //3D
	};
	//! @name Mesh functions
	//@{
	//! Meshes a set of line segments as constraints
	CDT_2 mesh2D(const std::vector<Mesh_Segment_2>& inSegVec, const MesherData& meshParams);
	//! Meshes a set of polygons as constraints
	CDT_2 mesh2D(const std::vector<std::vector<Mesh_Segment_2> >& inSegVec, const MesherData& meshParams);
	//! Meshes a set of polygons as constraints with material properties
	/*! @param inMatVec Contains either 1 or 2 GenericMaterial objects defining material properties for the entire
	 *  mesh or material properties for elements inside and outside the defining constraints
	 */
	CDT_2 mesh2D(const std::vector<std::vector<Mesh_Segment_2> >& inSegVec, const std::vector<GenericMaterial>& inMatVec, 
								const MesherData& meshParams);
	//! Meshes a set of line segments as constraints with material properties
	/*! @param inMatVec Contains either 1 or 2 GenericMaterial objects defining material properties for the entire
	 *  mesh or material properties for elements inside and outside the defining constraints
	 */
	CDT_2 mesh2D(const std::vector<Mesh_Segment_2>& inSegVec, const std::vector<GenericMaterial>& inMatVec,
								const MesherData& meshParams);
	//! Meshes a pixelated region with dimensions nx, ny and size width, height using triangular elements
	CDT_2 mesh2Dpixel(unsigned nx, unsigned ny,  double width, double height, const MesherData& meshParams);
	//! Meshes a pixelated region with dimensions nx, ny and size width, height using triangular elements and with material properties
	 /*!  @param inMatVec A vector of GenericMaterial objects, one for each pixel
	 */
	CDT_2 mesh2Dpixel( unsigned nx, unsigned ny,	double width, double height, const std::vector<GenericMaterial>& inMatVec, 
			const MesherData& meshParams);
	//! Meshes a pixelated region with dimensions nx, ny and size width, height using triangular elements and with material properties
	/*! @param inOptValVec A vector of optimization parameters, one for each pixel
	 */
	CDT_2 mesh2Dpixel(unsigned nx, unsigned ny, double width, double height, 
										const std::vector<double>& inOptValVec, const MesherData& meshParams);
	//! Meshes a region defined by the polyhedron @param 
	C3t3 mesh3D(const Mesh_polyhedron_3& polyToMesh, const MesherData& meshParams);
	//! Class to convert a TOMesh3DSurface object to a CGAL Polyhedron_3
	class PolyhedronBuilderFromTOMesh : public CGAL::Modifier_base<HalfedgeDS>
	{
	public:
		PolyhedronBuilderFromTOMesh(TOMesh3DSurface tom);
		void operator()( HalfedgeDS& hds);
	private:
		std::unique_ptr<TOMesh3DSurface> upTOM;
	};
	//@}
	//! @name Boundary mesh functions
	//@{
	//! Computes the signed area of the polygon defined by inSegVec
	Real computeSignedArea(const std::vector<Mesh_Segment_2>& inSegVec);
	//! Generates a set of line segments and inserts them into segVec
	void addPointsToSegments(const std::vector<Point_2_base>& ptVec, std::vector<Mesh_Segment_2>& segVec);
	//! Generates a set of line segments and inserts them into segVV as a new polygon
	void addPointsToSegments(const std::vector<Point_2_base>& ptVec, std::vector<std::vector<Mesh_Segment_2> >& segVV);
	//! Orders an arbitrary collection of line segments into a set of well-defined polygons (if possible)
	/*! Note: This uses an absolute tolerance. */
	void orderMeshSegments(const std::vector<Mesh_Segment_2>& segVec, std::vector<std::vector<Mesh_Segment_2> >& orderedVecVec);
	//! Copies the nested vector representation of line segments (collection of polygons) to one vector
	void flattenMeshSegments(const std::vector<std::vector<Mesh_Segment_2> >& orderedVec, std::vector<Mesh_Segment_2>& outSegVec);
	//! Check that mesh segments represent a valid boundary, i.e. no segments connecting at 3 or more points
	bool checkValidityAsBoundaryMesh(const std::vector<Mesh_Segment_2>& segVec);
	//! Creates a set of line segments from a list of points
	void generateSegment2FromPoints(const std::vector<Point_2_base>& ptVec, std::vector<Mesh_Segment_2>& segVec);
	//! Converts a set of line segments to a set of unique points
	void generatePointsFromSegment2(const std::vector<Mesh_Segment_2>& segVec, std::vector<Point_2_base>& ptVec);
	//! Converts a set of line segments into a CGAL polygon
	void meshSegs2Polygon(const std::vector<Mesh_Segment_2>& inSegVec, CGAL::Polygon_2<Mesh_K>& outPoly);
	//! Converts a set of line segments into several CGAL polygons: An outer boundary and any inner holes
	void meshSegs2PolygonWithHoles(const std::vector<Mesh_Segment_2>& inSegVec, CGAL::Polygon_2<Mesh_K>& outPoly, 
																std::vector<CGAL::Polygon_2<Mesh_K> >& holeVec);
	//! Converts a CGAL polygon into a vector of line segments
	void polygon2MeshSegs(const CGAL::Polygon_2<Mesh_K>& inPoly, std::vector<Mesh_Segment_2>& outSegVec);
	//! Projects chkPt to the boundary of the polygon defined in poly
	/*! First, the method finds the closest line segment of the polygon.  Next, an infinite supporting line
	 *  is formed and the orthogonal projection of the point onto the line is computed.  At this point, the projection
	 *  may still be outside of the polygon, so a check with tolerance tol is made, and if it fails, the projectsion
	 *  is determined to be the closest vertex of the polygon. 
	 */
	void projectPointToPolygon(Point_2_base& chkPt, const CGAL::Polygon_2<Mesh_K>& poly, double tol);
	//! Converts a polygon into a vector of coordinates.  
	/*! @param outComponentVec On return, this contains the points in the order: [p[0].x p[0].y ... p[N-1].x() p[N-1].y()]
	 *  assuming that the polygon's vertices are ordered as p[0] through p[N-1].
	 */
	void polygon2ComponentVec(const CGAL::Polygon_2<Mesh_K>& inPoly, std::vector<double>& outComponentVec);
	//! Converts a vector of coordinates to a polygon
	/*! @param inSegVec Vector of coordinates in the order: [p[0].x p[0].y ... p[N-1].x() p[N-1].y()]
	 */
	void componentVec2Polygon(const std::vector<double>& inSegVec, CGAL::Polygon_2<Mesh_K>& outPoly);
	//! Computes an approximate convex partition of inSegVec
	/*! @param outSegVec On return, contains a vector of vector of line segments, which define individual convex polygons.
	 */
	void partitionMeshSegPolys(const std::vector<Mesh_Segment_2>& inSegVec, std::vector<std::vector<Mesh_Segment_2> >& outSegVec);
	//@}
	//! @name 2D to 3D
	//@{
	//! Returns a linear extrusion of a 2D triangular mesh in TOMesh format
	TOMesh3DSurface linearExtrude(const CDT_2& inMesh, double length);
	//@}
}
}
#endif

