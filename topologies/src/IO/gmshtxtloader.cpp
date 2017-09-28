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

#include "gmshtxtloader.h"
#include "tomesh.h"
#include <iostream>
#include <string>
#include <sstream>
#include "cgal_types.h"

namespace Topologies{
namespace InputLoader
{
	namespace GMSHTxtMeshLoader
	{
		namespace
		{
			bool advanceStreamToPattern(const std::string& pattern, std::string& line, std::ifstream& inFile)
			{
				bool found = false;
				std::getline(inFile, line);
				while(!found && inFile.good())
				{
					if(line.find(pattern) != std::string::npos)
						found = true;
					else
						std::getline(inFile, line);
				}
				if(!found)
					std::cout << "ERROR: Couldn't parse exotxt file, pattern: " << pattern << " not found." << std::endl;
				return found;
			}

			bool isBoundaryElem(unsigned elemType, unsigned dim)
			{
				if(dim == 2 && (elemType == 1)) // 2D line
					return true;
				if(dim == 3 && (elemType == 2 || elemType == 3)) // 3d 
					return true;
				return false;
			}
			bool isSolidElem(unsigned elemType, unsigned dim)
			{
				if(dim == 2 && (elemType == 2 || elemType == 3)) // 2D tri or quad
					return true;
				if(dim == 3 && (elemType == 4 || elemType == 5)) // 3d tet or hex
					return true;
				return false;
			}
			unsigned getElementDimension(unsigned elemType)
			{
				// Only for low order elements
				if(elemType == 15) // Point
					return 0;
				else if(elemType == 1) // Line
					return 1;
				else if(elemType == 2 || elemType == 3) // Tri or quad
					return 2;
				return 3;
			}
			unsigned getNumNodesPerElement(unsigned elemType)
			{
				if(elemType == 3 || elemType == 4) // Quad or tet
					return 4;
				else if(elemType == 5) // Hex
					return 8;
				else if(elemType == 2) // Tri
					return 3;
				else if(elemType == 1) // Line
					return 2;
				else if(elemType == 15) // Point
					return 1;
				return 0;
			}
			std::pair<unsigned, std::vector<std::size_t>> read1Elem(std::ifstream& inFile, unsigned elemID)
			{
				// Get tags
				unsigned ntags;
				inFile >> ntags;
				unsigned pe = 1;
				if(ntags >= 1)
				{
					inFile >> pe;
					unsigned na;
					for(unsigned kt = 0; kt < ntags-1; ++kt)
						inFile >> na; // Not used
				}
				// Get connectivity
				std::vector<std::size_t> curCon;
				curCon.resize(getNumNodesPerElement(elemID));
				for(unsigned kn = 0; kn < curCon.size(); ++kn)
				{
					inFile >> curCon[kn];
					curCon[kn] -= 1; // Starts at 1
				}
				return std::make_pair(pe, curCon);
			}
			void readElements(std::ifstream& inFile, unsigned dim, std::vector<std::vector<std::size_t>>& elemVec, 
				std::vector<int>& tagIDVec)
			{
				std::size_t nelems;
				inFile >> nelems;
				elemVec.reserve(nelems);
				tagIDVec.reserve(nelems);
				for(std::size_t k = 0; k < nelems; ++k)
				{
					std::size_t id;
					inFile >> id; // Ignored
					unsigned elemID;
					inFile >> elemID;
					std::pair<unsigned, std::vector<std::size_t>> elemDef = read1Elem(inFile, elemID);
					if(isSolidElem(elemID, dim)) // Only save solid elements (boundaries are read in loadBoundaryPhysicalEntity)
					{
						elemVec.push_back(elemDef.second);
						tagIDVec.push_back(elemDef.first);
					}
				}
				elemVec.shrink_to_fit();
				tagIDVec.shrink_to_fit();
			}
		}

		std::unique_ptr<TOMesh3D> loadGMSHTxtFile3D(const std::string& fileName)
		{
			std::ifstream inFile(fileName.c_str());
			// Get coordinates
			std::string str;
			if(!advanceStreamToPattern("$Nodes", str, inFile))
				return nullptr;
			std::size_t nnodes;
			inFile >> nnodes;
			std::vector<Point_3_base> ptVec(nnodes);
			for(std::size_t k = 0; k < nnodes; ++k)
			{
				std::size_t id;
				inFile >> id; // Ignored, hopefully IDs are all sequential!
				double x, y, z;
				inFile >> x >> y >> z;
				ptVec[k] = Point_3_base(x, y, z);
			}
			// Read elements
			if(!advanceStreamToPattern("$Elements", str, inFile))
				return nullptr;
			std::vector<std::vector<std::size_t>> elemVec;
			std::vector<int> tagIDVec;
			readElements(inFile, 3, elemVec, tagIDVec);
			inFile.close();
			return std::unique_ptr<TOMesh3D>(new TOMesh3D(ptVec, elemVec, tagIDVec));
		}
		std::unique_ptr<TOMesh2D> loadGMSHTxtFile2D(const std::string& fileName)
		{
			std::ifstream inFile(fileName.c_str());
			// Get coordinates
			std::string str;
			if(!advanceStreamToPattern("$Nodes", str, inFile))
				return nullptr;
			std::size_t nnodes;
			inFile >> nnodes;
			std::vector<Point_2_base> ptVec(nnodes);
			for(std::size_t k = 0; k < nnodes; ++k)
			{
				std::size_t id;
				inFile >> id; // Ignored, hopefully IDs are all sequential!
				double x, y, z;
				inFile >> x >> y >> z;
				ptVec[k] = Point_2_base(x, y);
			}
			// Read elements
			if(!advanceStreamToPattern("$Elements", str, inFile))
				return nullptr;
			std::vector<std::vector<std::size_t>> elemVec;
			std::vector<int> tagIDVec;
			readElements(inFile, 2, elemVec, tagIDVec);
			inFile.close();
			return std::unique_ptr<TOMesh2D>(new TOMesh2D(ptVec, elemVec, tagIDVec));
		}
		std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> loadBoundaryPhysicalEntity(const std::string& fileName,
				unsigned dim)
		{
			std::ifstream inFile(fileName.c_str());
			std::string str;
			std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> boundaryPE;
			if(!advanceStreamToPattern("$Elements", str, inFile))
				return boundaryPE;
			std::size_t nelems;
			inFile >> nelems;
			// Read boundary elements into a map
			std::map<unsigned, std::vector<std::vector<std::size_t>>> boundaryMap;
			for(std::size_t k = 0; k < nelems; ++k)
			{
				std::size_t id;
				inFile >> id; // Ignored
				unsigned elemID;
				inFile >> elemID;
				std::pair<unsigned, std::vector<std::size_t>> elemDef = read1Elem(inFile, elemID);
				if(getElementDimension(elemID) < dim)
					boundaryMap[elemDef.first].push_back(elemDef.second);
			}
			// Store in expected format
			for(auto it = boundaryMap.begin(); it != boundaryMap.end(); ++it)
				boundaryPE.push_back(std::make_pair(it->first, it->second));
			return boundaryPE;
		}
	}
}
}
