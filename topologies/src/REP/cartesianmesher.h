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

#ifndef CARTESIANMESHER
#define CARTESIANMESHER

#include "cgal_types.h"
#include "tomesh.h"
#include "inputloader.h"
#include <vector>
#include <memory>

namespace Topologies{
//! Namespace for functions used in meshing Cartesian regions in 2d and 3d
namespace CartesianMesher
{
	//! Mesh a 2d rectangular region and return result in a TOMesh
	/*! MeshElementType can be either tri or quad, giving triangular elements or quadrilaterals. */
	std::unique_ptr<TOMesh2D> generateMesh(MeshElementType inMET, unsigned nx, unsigned ny, double width, double height, 
																			const std::vector<double>& optVals);
	//! Mesh a 3d rectangular solid region and return result in a TOMesh
	/*! MeshElementType can either be tet or hex, giving tetrahedral or hexahedral elements.*/
	std::unique_ptr<TOMesh3D> generateMesh(MeshElementType inMET, unsigned nx, unsigned ny, unsigned nz, 
																			double width, double length, double height,
                                      const std::vector<double>& optVals);
}
}

#endif

