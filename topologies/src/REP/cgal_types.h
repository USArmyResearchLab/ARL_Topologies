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

#ifndef CGAL_TYPES_H
#define CGAL_TYPES_H

// CGAL 4.10 fix
// To avoid error in CGAL/Mesh_criteria_3.h:
#define BOOST_PARAMETER_MAX_ARITY 12

/*! @file cgal_types.h
 *  A file containing typedefs and includes for CGAL
 */

// CGAL includes
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polygon_2_algorithms.h>

namespace Topologies{
// CGAL typedefs
typedef CGAL::Gmpq                                                  FT;
typedef CGAL::Extended_cartesian<FT>                                Extended_K;
typedef CGAL::Simple_cartesian<FT>                                  K;
typedef CGAL::Polygon_2<K>                                          Polygon_2;
typedef K::Segment_2                                                Segment_2;
typedef K::Vector_2                                                 Vector_2;
typedef K::Point_2                                                  Point_2;
typedef CGAL::Nef_polyhedron_2<Extended_K>                          Nef_polyhedron_2;
typedef Polygon_2::Vertex_iterator                                  Vertex_iterator;
typedef Polygon_2::Edge_const_circulator                            HE_circulator;
typedef Nef_polyhedron_2::Explorer                                  Explorer_2;
}

// CGAL includes 3D
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Gmpq.h>

namespace Topologies{
// CGAL typedefs
typedef CGAL::Convex_hull_traits_3<K>                               CH_traits_3;
typedef CH_traits_3::Polyhedron_3                                   Polyhedron_3;
typedef CGAL::Plane_3<K>                                            Plane_3;
typedef K::Segment_3                                                Segment_3;
typedef K::Vector_3                                                 Vector_3;
typedef K::Point_3                                                  Point_3;
typedef CGAL::Nef_polyhedron_3<K>                                   Nef_polyhedron_3;
typedef Polyhedron_3::Vertex_iterator                               Vertex_iterator_3;
typedef Polyhedron_3::Facet_iterator                                Facet_iterator_3;
typedef Polyhedron_3::Point_iterator                                Point_iterator_3;
typedef Polyhedron_3::Facet::Halfedge_around_facet_const_circulator HE_circulator_3;
}

// 2D Meshing
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include "genericmaterial.h"

namespace Topologies{
//! A struct for use in CGAL's mesher.  This allows placing more information in each mesh face.
struct FaceInfo2
{
  FaceInfo2() : nesting_level(0), optVal(0.) {}
  int nesting_level;
  GenericMaterial triMat;
  double optVal;

  bool in_domain(){
    return nesting_level%2 == 1;
  }
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel             Mesh_K;
typedef CGAL::Triangulation_vertex_base_with_info_2<int,Mesh_K>         Vbb;
typedef CGAL::Triangulation_vertex_base_2<Mesh_K,Vbb>                   Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Mesh_K>     Fbbb;
typedef CGAL::Constrained_triangulation_face_base_2<Mesh_K, Fbbb>       Fbb;
typedef CGAL::Delaunay_mesh_face_base_2<Mesh_K,Fbb>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                    Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<Mesh_K, Tds>         CDT_2;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT_2>                      Criteria;
typedef CGAL::Delaunay_mesher_2<CDT_2, Criteria>                        Mesher_2;
}

// 2D alpha shapes
#include <CGAL/Weighted_point.h>
#include <CGAL/Weighted_alpha_shape_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

namespace Topologies{
typedef Mesh_K::FT                            Mesh_FT;
typedef Mesh_K::Point_2                         Point_2_base;
typedef Mesh_K::Point_3							Point_3_base;
typedef CGAL::Weighted_point<Point_2_base, Mesh_FT>           W_Point_2;
typedef Mesh_K::Segment_2                       Mesh_Segment_2;
typedef Mesh_K::Vector_2                        Mesh_Vector_2;

typedef CGAL::Weighted_alpha_shape_euclidean_traits_2<Mesh_K>     Gt;
typedef CGAL::Regular_triangulation_vertex_base_2<Gt>           Rvb;
typedef CGAL::Alpha_shape_vertex_base_2<Gt,Rvb>             aVb;
typedef CGAL::Regular_triangulation_face_base_2<Gt>           Rf;
typedef CGAL::Alpha_shape_face_base_2<Gt, Rf>               aFb;

typedef CGAL::Triangulation_data_structure_2<aVb,aFb>           aTds;
typedef CGAL::Regular_triangulation_2<Gt,aTds>              Triangulation_2;

typedef CGAL::Alpha_shape_2<Triangulation_2>                Alpha_shape_2;

typedef Alpha_shape_2::Face                       Face;
typedef Alpha_shape_2::Vertex                       Vertex;
typedef Alpha_shape_2::Edge                       Edge;
typedef Alpha_shape_2::Face_handle                    Face_handle;
typedef Alpha_shape_2::Vertex_handle                  Vertex_handle;

typedef Alpha_shape_2::Face_circulator                  Face_circulator;
typedef Alpha_shape_2::Vertex_circulator                  Vertex_circulator;

typedef Alpha_shape_2::Locate_type                    Locate_type;

typedef Alpha_shape_2::Face_iterator                    Alpha_Face_iterator;
typedef Alpha_shape_2::Vertex_iterator                  Alpha_Vertex_iterator;
typedef Alpha_shape_2::Edge_iterator                    Alpha_Edge_iterator;
typedef Alpha_shape_2::Edge_circulator                  Alpha_Edge_circulator;

typedef Alpha_shape_2::Alpha_iterator                   Alpha_iterator;
typedef Alpha_shape_2::Alpha_shape_edges_iterator             Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator                    Alpha_shape_vertices_iterator;
}

// 3D meshing

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
//#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace Topologies{
// Domain
//typedef CGAL::Polyhedron_3<Mesh_K>										Mesh_polyhedron_3;
//typedef Mesh_polyhedron_3::HalfedgeDS      								HalfedgeDS;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Mesh_K>			Mesh_domain;
typedef CGAL::Mesh_polyhedron_3<Mesh_K>::type									Mesh_polyhedron_3;	
typedef Mesh_polyhedron_3::HalfedgeDS									HalfedgeDS;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type 					Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, 
		Mesh_domain::Corner_index, Mesh_domain::Curve_segment_index> 	C3t3;
typedef CGAL::Mesh_criteria_3<Tr> 										Mesh_criteria;
}
// 3D Surface meshing
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
namespace Topologies{
typedef CGAL::Surface_mesh_default_triangulation_3 SM3_Tr;
typedef CGAL::Complex_2_in_triangulation_3<SM3_Tr> C2t3;
typedef SM3_Tr::Geom_traits Tr_GT;
}
#endif

