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

#include "stltxtloader.h"
#include "cgal_types.h"
#include "tomesh.h"
#include "filter3d.h" // For Point_3_hash
#include "geometrytranslation.h"
#include <vector>
#include <unordered_map>
#include <iostream>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace Topologies{
namespace InputLoader
{
	namespace STLTxtLoader
	{
		namespace
		{
			bool readSTL(std::ifstream& inFile, std::vector<Point_3_base>& ptVec, std::vector<std::vector<std::size_t>>& connectivityVec)
			{
				// Define stl format strings
				std::string solidStr("solid"), facetStr("facet"), outerStr("outer"), loopStr("loop"), vertexStr("vertex"),
					normalStr("normal"), endLoopStr("endloop"), endFacetStr("endfacet"), endsolidStr("endsolid");
				std::string curToken;
				// Data structures
				std::unordered_map<Point_3_base, std::size_t, Point_3_hash> ptIDMap; // Using a map to define connectivity
				// Start reading
				inFile >> curToken;
				if(curToken != solidStr)
					return false;
				std::getline(inFile, curToken); // Read and discard solid name
				inFile >> curToken; // Should be "facet normal" or "endsolid"
				while(curToken != endsolidStr)
				{
					if(curToken != facetStr)
						return false;
					std::getline(inFile, curToken); // Read and discard normal data
					inFile >> curToken; // Should be "outer"
					if(curToken != outerStr)
						return false;
					inFile >> curToken; // Should be "loop"
					if(curToken != loopStr)
						return false;
					std::vector<std::size_t> newTri(3);
					for(unsigned k = 0; k < 3; ++k)
					{
						inFile >> curToken; // Should be "vertex"
						if(curToken != vertexStr)
							return false;
						Point_3_base p;
						inFile >> p;
						auto insertionResult = ptIDMap.insert(std::make_pair(p, 0)); // Second value is irrelevant, will be set below
						if(insertionResult.second) // New insertion
						{
							ptVec.push_back(p); // Store point
							ptIDMap[p] = ptVec.size() - 1; // Assign correct id and increment
						}
						newTri[k] = ptIDMap[p];
					}
					inFile >> curToken; // Should be "endloop"
					if(curToken != endLoopStr)
						return false;
					inFile >> curToken; // Should be "endfacet"
					if(curToken != endFacetStr)
						return false;
					connectivityVec.push_back(newTri);
					inFile >> curToken; // Should be "facet normal" or "endsolid"
				}
				return true; // No error
			}
		}
		
		Mesh_polyhedron_3 loadFileCGAL(const std::string& fileName)
		{
			// Read STL into TOMesh format
			std::unique_ptr<TOMesh3DSurface> upTOM = loadFileTOMesh(fileName);
			if(!upTOM)
				return Mesh_polyhedron_3(); // Error reading STL file
			// Convert to CGAL polyhedron using incremental builder
			Mesh_polyhedron_3 MP3;
			GeometryTranslation::PolyhedronBuilderFromTOMesh builder(*upTOM);
			MP3.delegate(builder);
			return MP3;
		}

		std::unique_ptr<TOMesh3DSurface> loadFileTOMesh(const std::string& fileName)
		{
			std::ifstream inFile(fileName.c_str());
			std::vector<Point_3_base> ptVec;
			std::vector<std::vector<std::size_t>> connectivityVec;
			bool success = readSTL(inFile, ptVec, connectivityVec);
			if(!success)
				return nullptr;
			return std::unique_ptr<TOMesh3DSurface>(new TOMesh3DSurface(std::move(ptVec), std::move(connectivityVec)));
		}
	}
}
}

