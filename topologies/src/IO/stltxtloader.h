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

#ifndef STLTXTINPUTLOADER_H
#define STLTXTINPUTLOADER_H

#include "inputloader.h"
#include "cgal_types.h"
#include <memory>

namespace Topologies{
class TOMesh3DSurface;

namespace InputLoader
{
	//! A collection of functions to load ASCII STL files 
	namespace STLTxtLoader
	{
		//! Loads a 3d ExodusII files and returns the mesh in a TOMesh3D object
		Mesh_polyhedron_3 loadFileCGAL(const std::string& fileName);
		//! Loads a 3d ExodusII files and returns the mesh in a TOMesh2D object
		std::unique_ptr<TOMesh3DSurface> loadFileTOMesh(const std::string& fileName);
	}
}
}
#endif

