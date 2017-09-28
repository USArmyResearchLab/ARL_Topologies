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

#ifndef OUTPUTWRITER_H
#define OUTPUTWRITER_H

#include <string>
#include <vector>
#include "cgal_types.h"
#include "postprocess.h"

namespace Topologies{
class TOMesh;
class TopOptRep;

//! A collection of functions to facilitate the outputing of meshes in various formats
namespace OutputWriter
{
	//! @name Meshes
	//@{
	//! Write @param inMesh to stream @param outFile in Matlab format
	void plotMeshMatlab(const TOMesh* const inMesh, std::ofstream& outFile);
	//! Write @param inMesh to the file named @param fileName in Matlab format
	void plotMeshMatlab(const TOMesh* const inMesh, const std::string& fileName);
	//! Write @param inMesh to file @param fileName in STL format in ascii or binary (@param useASCII)
	void plotSTL(const TOMesh* const inMesh, const std::string& fileName, bool useASCII = true);
	//! Write @param inMesh to file @param fileName in STL format in ascii or binary (@param useASCII)
	void plotSTL(const TOMesh& inMesh, const std::string& fileName, bool useASCII = true);
	//! Write @param surfMesh to file @param fileName in STL format in ascii or binary (@param useASCII)
	void plotSTL(const C2t3& surfMesh, const std::string& fileName, bool useASCII = true);
	//! Write line segment list @param setVV to file @param fileName in Matlab format
	void plotSegmentsMatlab(const std::vector<std::vector<Mesh_Segment_2> >& segVV, const std::string& fileName);
	//! Write line segment list @param segVec to file @param fileName in Matlab format
	void plotSegmentsMatlab(const std::vector<Mesh_Segment_2>& segVec, const std::string& fileName);
	//! Write @paraminMesh to stream @param outFile in VTK format
	void plotMeshVTK(const TOMesh* const inMesh, std::ofstream& outFile);
	//! Write @param inMesh to file @param fileName in VTK format
	void plotMeshVTK(const TOMesh* const inMesh, const std::string& fileName);
	//! Write @param inMesh to file @param fileName in VTK format and append number `kt` to file name
	void plotMeshVTK(const TOMesh* const inMesh, unsigned kt, const std::string& fileName);
	//! Write @param inMesh to stream @param outFile in GMSH format.
	//! @param writeMesh Write entire mesh on true, just data on false
	void plotMeshGMSH(const TOMesh* const inMesh, std::ofstream& outFile, unsigned kt = 0, bool writeMesh = true);
	//! Write @param inMesh to file @param fileName in GMSH format and label time step as `kt`
	//! @param writeMesh Write entire mesh on true, just data on false
	void plotMeshGMSH(const TOMesh* const inMesh, const std::string& fileName, unsigned kt = 0, bool writeMesh = true);
	//@}
	//! @name Data
	//@{
	//! Output raw data in Matlab format for a 3d array
	void plotDataMatlab(const std::vector<double>& data, std::size_t nx, std::size_t ny, std::size_t nz, const std::string& fileName);
	//! Output raw data in Matlab format for a 2d array
	void plotDataMatlab(const std::vector<double>& data, std::size_t nx, std::size_t ny, const std::string& fileName);
	//! Output raw data of `outTOR` to file `fileName`
	void outputTORData(const TopOptRep* outTOR, const std::string& fileName);
	//@}
}
}
#endif

