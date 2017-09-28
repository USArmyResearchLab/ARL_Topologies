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

#include "exotxtmeshloader.h"
#include "tomesh.h"
#include <iostream>
#include <string>
#include <sstream>
#include "cgal_types.h"

namespace Topologies{
namespace InputLoader
{
	namespace ExoTxtMeshLoader
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
			
			std::vector<unsigned> getTetFaceVertices(unsigned faceID)
			{
				// Hard-coded vertex ids, as observed from the exodus file format
				std::vector<unsigned> vertexIDs;
				if(faceID == 1)
					vertexIDs = {1, 2, 4};
				else if(faceID == 2)
					vertexIDs = {2, 3, 4};
				else if(faceID == 3)
					vertexIDs = {1, 4, 3};
				else if(faceID == 4)
					vertexIDs = {1, 3, 2};
				return vertexIDs;
			}

			std::vector<unsigned> getHexFaceVertices(unsigned faceID)
			{
				// Hard-coded vertex ids, as observed from the exodus file format
				std::vector<unsigned> vertexIDs;
				if(faceID == 1)
					vertexIDs = {2, 1, 5, 6};
				else if(faceID == 2)
					vertexIDs = {2, 3, 7, 6};
				else if(faceID == 3)
					vertexIDs = {3, 4, 8, 7};
				else if(faceID == 4)
					vertexIDs = {1, 5, 8, 4};
				else if(faceID == 5)
					vertexIDs = {2, 1, 4, 3};
				else if(faceID == 6)
					vertexIDs = {5, 6, 7, 8};
				return vertexIDs;
			}

			std::vector<unsigned> getTriEdgeVertices(unsigned edgeID)
			{
				// Hard-coded vertex ids, as observed from the exodus file format
				std::vector<unsigned> vertexIDs;
				if(edgeID == 1)
					vertexIDs = {1, 2};
				else if(edgeID == 2)
					vertexIDs = {2, 3};
				else if(edgeID == 3)
					vertexIDs = {3, 1};
				return vertexIDs;
			}

			std::vector<unsigned> getQuadEdgeVertices(unsigned edgeID)
			{
				// Hard-coded vertex ids, as observed from the exodus file format
				std::vector<unsigned> vertexIDs;
				if(edgeID == 1)
					vertexIDs = {1, 2};
				else if(edgeID == 2)
					vertexIDs = {2, 3};
				else if(edgeID == 3)
					vertexIDs = {3, 4};
				else if(edgeID == 4)
					vertexIDs = {4, 1};
				return vertexIDs;
			}

			unsigned getNumNodesPerElement(const std::string& elemType)
			{
				if(elemType == "TETRA" || elemType == "SHELL4" || elemType == "QUAD4")
					return 4;
				else if(elemType == "HEX8")
					return 8;
				else if(elemType == "TRI3")
					return 3;
				return 0;
			}
		}

		unsigned readMeshDimension(const std::string& fileName)
		{
			std::ifstream inFile(fileName.c_str());
			std::string str;
			if(!advanceStreamToPattern("! dimensions", str, inFile))
				return 0;
			std::stringstream ss1(str);
			unsigned ndim;
			ss1 >> ndim;
			inFile.close();
			return ndim;
		}

		std::unique_ptr<TOMesh3D> loadExoTxtFileTet(const std::string& fileName)
		{
			std::ifstream inFile(fileName.c_str());
			// Get coordinates
			std::string str;
			if(!advanceStreamToPattern("! dimensions", str, inFile))
				return nullptr;
			std::stringstream ss1(str);
			unsigned ndim;
			ss1 >> ndim;
			if(ndim != 3)
			{
				std::cout << "Error reading file: " << fileName << ", expected 3 dimensions, got " << ndim << std::endl;
				return nullptr;
			}
			if(!advanceStreamToPattern("! nodes, elements", str, inFile))
				return nullptr;
			// Have line with number of nodes and elements in str
			std::stringstream ss2(str);
			std::size_t nnodes, nelems, nblocks;
			ss2 >> nnodes >> nelems >> nblocks;
			// Read lines until we get to the coordinates
			if(!advanceStreamToPattern("! Coordinates", str, inFile))
				return nullptr;
			// Read coords
			std::vector<Point_3_base> ptVec(nnodes);
			for(std::size_t k = 0; k < nnodes; ++k)
			{
				double x, y, z;
				inFile >> x >> y >> z;
				ptVec[k] = Point_3_base(x, y, z);
			}
			// Read connectivity
			std::vector<std::vector<std::size_t>> connectivityVec(nelems);
			std::vector<int> blockIDVec(nelems);
			std::vector<int>::iterator itBID = blockIDVec.begin();
			std::size_t kelem = 0;
			for(std::size_t kb = 0; kb < nblocks; ++kb)
			{
				if(!advanceStreamToPattern("! Element block", str, inFile))
					return nullptr;
				if(!advanceStreamToPattern("! ID, elements, element type", str, inFile))
					return nullptr;
				unsigned blockID;
				std::size_t nelemsblock;
				std::string elemType;
				std::stringstream ss3(str);
				ss3 >> blockID >> nelemsblock >> elemType;
				// Assign block ids
				std::fill(itBID, itBID + nelemsblock, blockID);
				itBID += nelemsblock;
				unsigned nnodesperelem = getNumNodesPerElement(elemType);
				if(nnodesperelem == 0)
				{
					std::cout << "Unknown element type found in exodus file: " << elemType << std::endl;
					return nullptr;
				}
				if(!advanceStreamToPattern("! Connectivity", str, inFile))
					return nullptr;
				// Read elements
				for(std::size_t ke = 0; ke < nelemsblock; ++ke)
				{
					std::vector<std::size_t>& curCon = connectivityVec[kelem++];
					curCon.resize(nnodesperelem);
					for(unsigned kn = 0; kn < nnodesperelem; ++kn)
					{
						inFile >> curCon[kn];
						curCon[kn] -= 1;
					}
				}
			}
			// Set up TOMesh
			inFile.close();
			return std::unique_ptr<TOMesh3D>(new TOMesh3D(ptVec, connectivityVec, blockIDVec));
		}

		std::unique_ptr<TOMesh2D> loadExoTxtFileTri(const std::string& fileName)
		{
			std::ifstream inFile(fileName.c_str());
			// Get coordinates
			std::string str;
			if(!advanceStreamToPattern("! dimensions", str, inFile))
				return nullptr;
			std::stringstream ss1(str);
			unsigned ndim;
			ss1 >> ndim;
			if(ndim != 2)
			{
				std::cout << "Error reading file: " << fileName << ", expected 2 dimensions, got " << ndim << std::endl;
				return nullptr;
			}
			if(!advanceStreamToPattern("! nodes, elements", str, inFile))
				return nullptr;
			// Have line with number of nodes and elements in str
			std::stringstream ss2(str);
			std::size_t nnodes, nelems, nblocks;
			ss2 >> nnodes >> nelems >> nblocks;
			// Read lines until we get to the coordinates
			if(!advanceStreamToPattern("! Coordinates", str, inFile))
				return nullptr;
			// Read coords
			std::vector<Point_2_base> ptVec(nnodes);
			for(std::size_t k = 0; k < nnodes; ++k)
			{
        double x, y;
        inFile >> x >> y;
				ptVec[k] = Point_2_base(x, y);
			}
			// Read connectivity
			std::vector<std::vector<std::size_t>> connectivityVec(nelems);
			std::vector<int> blockIDVec(nelems);
			std::vector<int>::iterator itBID = blockIDVec.begin();
			std::size_t kelem = 0;
			for(unsigned kb = 0; kb < nblocks; ++kb)
			{
				if(!advanceStreamToPattern("! Element block", str, inFile))
					return nullptr;
				if(!advanceStreamToPattern("! ID, elements, element type", str, inFile))
					return nullptr;
				unsigned blockID;
				std::size_t nelemsblock;
				std::string elemType;
				std::stringstream ss3(str);
				ss3 >> blockID >> nelemsblock >> elemType;
				// Assign block ids
				std::fill(itBID, itBID + nelemsblock, blockID);
				itBID += nelemsblock;
				unsigned nnodesperelem = getNumNodesPerElement(elemType);
				if(nnodesperelem == 0) // Check for unknown element type
					return nullptr;
				if(!advanceStreamToPattern("! Connectivity", str, inFile))
					return nullptr;
				for(std::size_t ke = 0; ke < nelemsblock; ++ke)
				{
					std::vector<std::size_t>& curCon = connectivityVec[kelem++];
					curCon.resize(nnodesperelem);
					for(unsigned kn = 0; kn < nnodesperelem; ++kn)
					{
						inFile >> curCon[kn];
						curCon[kn] -= 1; // Starts at 1
					}
				}
			}
			// Set up TOMesh
			inFile.close();
			return std::unique_ptr<TOMesh2D>(new TOMesh2D(ptVec, connectivityVec, blockIDVec));
		}

		std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> loadSideSets(const std::string& fileName, 
																																											const TOMesh* const mesh)
		{
			std::vector<std::pair<unsigned, std::vector<std::vector<std::size_t>>>> sideSetVV;
			std::ifstream inFile(fileName.c_str());
			std::string str;
			if(!advanceStreamToPattern("! dimensions", str, inFile))
				return sideSetVV;
			std::stringstream ss1(str);
			unsigned ndim;
			ss1 >> ndim;
			if(!advanceStreamToPattern("! #node sets, #side sets", str, inFile))
				return sideSetVV;
			unsigned nnodesets, nsidesets;
			std::stringstream ss2(str);
			ss2 >> nnodesets >> nsidesets;
			sideSetVV.resize(nsidesets);
			if(!advanceStreamToPattern("! Number of Side Sets", str, inFile))
				return sideSetVV;
			for(unsigned k = 0; k < nsidesets; ++k)
			{
				if(!advanceStreamToPattern("! ID, #sides, #nodes", str, inFile))
					return sideSetVV;
				std::stringstream ss2(str);
				unsigned sideSetID;
				std::size_t nelems;
				ss2 >> sideSetID >> nelems;
				if(!advanceStreamToPattern("! Elements and sides for side sets", str, inFile))
					return sideSetVV;
				std::vector<std::vector<std::size_t>> curFaceVec(nelems);
				for(std::size_t ke = 0; ke < nelems; ++ke)
				{
					// Read element ID and face/line number
					unsigned elemID, faceID;
					inFile >> elemID >> faceID;
					elemID--;
					// Get corresponding vertex ids
					std::vector<std::size_t> curElem = mesh->getElementConnectivity(elemID);
					std::vector<unsigned> vertexIDVec;
					if(ndim == 2 && curElem.size() == 3)
						vertexIDVec = getTriEdgeVertices(faceID);
					else if(ndim == 2 && curElem.size() == 4)
						vertexIDVec = getQuadEdgeVertices(faceID);
					else if(ndim == 3 && curElem.size() == 4)
						vertexIDVec = getTetFaceVertices(faceID);
					else if(ndim == 3 && curElem.size() == 8)
						vertexIDVec = getHexFaceVertices(faceID);
					curFaceVec[ke].resize(vertexIDVec.size());
					for(std::size_t kv = 0; kv < vertexIDVec.size(); ++kv)
						curFaceVec[ke][kv] = curElem[vertexIDVec[kv] - 1];
				}
				sideSetVV[k] = std::make_pair(sideSetID, curFaceVec);
			}
			inFile.close();
			return sideSetVV;
		}

		std::vector<std::pair<unsigned, std::vector<std::size_t>>> loadNodeSets(const std::string& fileName)
		{
			std::vector<std::pair<unsigned, std::vector<std::size_t>>> nodeSetVV;
			std::ifstream inFile(fileName.c_str());
			std::string str;
			if(!advanceStreamToPattern("! #node sets, #side sets", str, inFile))
				return nodeSetVV;
			unsigned nnodesets;
			std::stringstream ss1(str);
			ss1 >> nnodesets;
			nodeSetVV.resize(nnodesets);
			for(unsigned k = 0; k < nnodesets; ++k)
			{
				if(!advanceStreamToPattern("! Nodal point set", str, inFile))
					return nodeSetVV;
				if(!advanceStreamToPattern("! ID, nodes", str, inFile))
					return nodeSetVV;
				std::stringstream ss2(str);
				unsigned nodeSetID;
				std::size_t nnodes;
				ss2 >> nodeSetID >> nnodes;
				std::vector<std::size_t> curNodeVec(nnodes);
				for(unsigned kn = 0; kn < nnodes; ++kn)
				{
					inFile >> curNodeVec[kn];
					curNodeVec[kn] -= 1; // Exodus orders from 1, TOMesh orders from 0
					inFile.ignore(512, '\n');
				}
				nodeSetVV[k] = std::make_pair(nodeSetID, curNodeVec);
			}
			inFile.close();
			return nodeSetVV;
		}

	}
}
}
