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

#ifndef EXOTXTINPUTLOADER_H
#define EXOTXTINPUTLOADER_H

#include "inputloader.h"
#include <memory>

namespace Topologies{
class TOMesh3D;
class TOMesh2D;
class TOMesh;

namespace InputLoader
{
	//! A collection of functions to load text-based Exodus II files
	/*! The Exodus II file format can be generated from Cubit and other pre-processing software.
	 *  These files only parse text-based Exodus II files, not binary.  The Trilinos library
	 *  has a tool to convert binary Exodus files to text called exotxt.  
	 */
	namespace ExoTxtMeshLoader
	{
		//! Loads a 3d ExodusII files and returns the mesh in a TOMesh3D object
		std::unique_ptr<TOMesh3D> loadExoTxtFileTet(const std::string& fileName);
		//! Loads a 2d ExodusII files and returns the mesh in a TOMesh2D object
		std::unique_ptr<TOMesh2D> loadExoTxtFileTri(const std::string& fileName);
		//! Returns the dimension of the mesh in a given file (2 or 3)
		unsigned readMeshDimension(const std::string& fileName);
		//! Loads 'side sets' (boundary elements) from an ExodusII files
		/*! A side set is the way an ExodusII file designates collections of boundary elements.
		 *  These are used to define boundary conditions.  The format use here to store side sets
		 *  is a vector of pairs that contain the side set ID and the boundary element connectivity.
		 *  The boundary element connectivity is a vector of vectors of unsigned ints, indicating node
		 *  IDs for each element.  
		 */
		std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> loadSideSets(const std::string& fileName, 
																							const TOMesh* const mesh);
		//! Loads ExodusII 'node sets', which are collections of nodes to specify boundary conditions.
		/*! Similar to a 'side set', a node set is a collection of nodes.  The format is a vector of pairs
		 *  of the node set ID and a vector of node IDs.  
		 */
		std::vector<std::pair<unsigned, std::vector<std::size_t>>> loadNodeSets(const std::string& fileName);
	}
}
}
#endif

