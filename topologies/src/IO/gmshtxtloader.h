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

#ifndef GMSHTXTINPUTLOADER_H
#define GMSHTXTINPUTLOADER_H

#include "inputloader.h"
#include <memory>

namespace Topologies{
class TOMesh3D;
class TOMesh2D;
class TOMesh;

namespace InputLoader
{
	//! A collection of functions to load text-based GMSH files
	/*! The msh file format can be generated from GMSH and other pre-processing software.
	 *  These files only parse text-based GMSH files, not binary.
	 */
	namespace GMSHTxtMeshLoader
	{
		//! Loads a 3d GMSH file and returns the mesh in a TOMesh3D object
		std::unique_ptr<TOMesh3D> loadGMSHTxtFile3D(const std::string& fileName);
		//! Loads a 2d GMSH file and returns the mesh in a TOMesh2D object
		std::unique_ptr<TOMesh2D> loadGMSHTxtFile2D(const std::string& fileName);
		//! Loads boundary element physical entity from GMSH file
		/*! A physical entity is the way GMSH designates collections of boundary elements.
		 *  These are used to define boundary conditions.  The format use here to store PEs
		 *  is a vector of pairs that contain the PE ID and the boundary element connectivity.
		 *  The boundary element connectivity is a vector of vectors of unsigned ints, indicating node
		 *  IDs for each element.
		 */
		std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> loadBoundaryPhysicalEntity(const std::string& fileName,
				unsigned dim);
	}
}
}
#endif

