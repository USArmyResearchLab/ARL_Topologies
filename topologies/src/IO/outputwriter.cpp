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

#include "outputwriter.h"
#include "cgal_types.h"
#include "tomesh.h"
#include "topoptrep.h"
#include <string>
#include <fstream>
#include <vector>
#include <CGAL/IO/output_surface_facets_to_triangle_soup.h>

namespace Topologies{
namespace OutputWriter
{
	namespace
  {
		template<typename numberType>
		numberType swapEndianness(numberType num)
		{
			numberType res = 0;
			std::size_t nbytes = sizeof(num);
			unsigned char* numBuf = (unsigned char*)&num;
			unsigned char* resBuf = (unsigned char*)&res;
			for(std::size_t k = 0; k < nbytes; ++k)
				resBuf[nbytes - k - 1] = numBuf[k];
			return res;
		}

		template<typename floatType, typename CGALType_3>
		void writePointBinary(CGALType_3 pt, std::ofstream& plotFile)
		{
			for(unsigned short k = 0; k < 3; ++k)
			{
				floatType coord = (floatType)pt[k];
				plotFile.write((char*)&coord, sizeof(coord));
			}
		}

		void plotSTLASCII(const TOMesh& inMesh, const std::string& fileName)
		{
			std::ofstream outFile(fileName);
			outFile.precision(8);
			outFile << std::scientific;
			std::cout << "Outputting STL, num triangles: " << inMesh.getNumElements() << '\n';
			outFile << "solid " << fileName << '\n';
			for(std::size_t ke = 0; ke < inMesh.getNumElements(); ++ke)
			{
				const std::vector<std::size_t> curElem = inMesh.getElementConnectivity(ke);
				CGAL::Point_3<Mesh_K> p0 = inMesh.getNode3D(curElem[0]), p1 = inMesh.getNode3D(curElem[1]),
					p2 = inMesh.getNode3D(curElem[2]);
				CGAL::Vector_3<Mesh_K> n = CGAL::normal(p0, p1, p2);
				outFile << "facet normal " << n.x() << " " << n.y() << " " << n.z() << "\n";
				outFile << " outer loop\n";
				outFile << "  vertex " << p0.x() << " " << p0.y() << " " << p0.z() << "\n";
				outFile << "  vertex " << p1.x() << " " << p1.y() << " " << p1.z() << "\n";
				outFile << "  vertex " << p2.x() << " " << p2.y() << " " << p2.z() << "\n";
				outFile << " endloop\nendfacet" << '\n';;
			}
			outFile << "endsolid" << '\n';
			outFile.close();
		}

		void plotSTLASCII(const C2t3& surfMesh, const std::string& fileName)
		{
			typedef std::vector<CGAL::Triangle_3<Mesh_K> > Triangle_vector;
			Triangle_vector triVec;
			CGAL::output_surface_facets_to_triangle_soup(surfMesh, std::back_inserter(triVec));
			std::ofstream outFile(fileName);
			outFile.precision(8);
			outFile << std::scientific;
			outFile << "solid " << fileName << '\n';
			for(Triangle_vector::iterator it = triVec.begin(); it != triVec.end(); ++it)
			{
				CGAL::Point_3<Mesh_K> p0 = it->vertex(0), p1 = it->vertex(1), p2 = it->vertex(2);
				CGAL::Vector_3<Mesh_K> n = CGAL::normal(p0, p1, p2);
				outFile << "facet normal " << n.x() << " " << n.y() << " " << n.z() << "\n";
				outFile << " outer loop\n";
				outFile << "  vertex " << fabs(p0.x()) << " " << fabs(p0.y()) << " " << fabs(p0.z()) << "\n";
				outFile << "  vertex " << fabs(p1.x()) << " " << fabs(p1.y()) << " " << fabs(p1.z()) << "\n";
				outFile << "  vertex " << fabs(p2.x()) << " " << fabs(p2.y()) << " " << fabs(p2.z()) << "\n";
				outFile << " endloop\nendfacet" << '\n';
			}
			outFile << "endsolid" << '\n';
			outFile.close();
		}

		void plotSTLBinary(const TOMesh& inMesh, const std::string& fileName)
		{
			std::ofstream outFile(fileName, std::ios::in|std::ios::binary|std::ios::trunc);
			assert(outFile.is_open());
			std::cout << "Outputting STL, num triangles: " << inMesh.getNumElements() << '\n';
			char header[80] = {0}; // 80 byte header, meaningless?
			outFile.write(header, 80);
			// Write number of facets
			uint32_t nfacets = swapEndianness((uint32_t)inMesh.getNumElements());
			outFile.write((char*)&nfacets, sizeof(nfacets)); // Should be 4 bytes
			for(std::size_t ke = 0; ke < inMesh.getNumElements(); ++ke)
			{
				const std::vector<std::size_t> curElem = inMesh.getElementConnectivity(ke);
				CGAL::Point_3<Mesh_K> p0 = inMesh.getNode3D(curElem[0]), p1 = inMesh.getNode3D(curElem[1]),
					p2 = inMesh.getNode3D(curElem[2]);
				CGAL::Vector_3<Mesh_K> n = CGAL::normal(p0, p1, p2);
				writePointBinary<float, CGAL::Vector_3<Mesh_K>>(n, outFile);
				writePointBinary<float, CGAL::Point_3<Mesh_K>>(p0, outFile);
				writePointBinary<float, CGAL::Point_3<Mesh_K>>(p1, outFile);
				writePointBinary<float, CGAL::Point_3<Mesh_K>>(p2, outFile);
				// Write attribute
				char buf[2] = {0, 0};
				outFile.write(buf, 2); // Set to zero
			}
			outFile.close();
		}

		void plotSTLBinary(const C2t3& surfMesh, const std::string& fileName)
		{
			std::ofstream outFile(fileName, std::ios::in|std::ios::binary|std::ios::trunc);
			assert(outFile.is_open());
			char header[80] = {0}; // 80 byte header, meaningless?
			outFile.write(header, 80);
			// Write number of facets
			typedef std::vector<CGAL::Triangle_3<Mesh_K>> Triangle_vector;
			Triangle_vector triVec;
			CGAL::output_surface_facets_to_triangle_soup(surfMesh, std::back_inserter(triVec));
			uint32_t nfacets = swapEndianness((uint32_t)triVec.size());
			outFile.write((char*)&nfacets, sizeof(nfacets)); // Should be 4 bytes
			for(Triangle_vector::const_iterator it = triVec.begin(); it != triVec.end(); ++it)
			{
				CGAL::Point_3<Mesh_K> p0 = it->vertex(0), p1 = it->vertex(1), p2 = it->vertex(2);
				CGAL::Vector_3<Mesh_K> n = CGAL::normal(p0, p1, p2);
				writePointBinary<float, CGAL::Vector_3<Mesh_K>>(n, outFile);
				writePointBinary<float, CGAL::Point_3<Mesh_K>>(p0, outFile);
				writePointBinary<float, CGAL::Point_3<Mesh_K>>(p1, outFile);
				writePointBinary<float, CGAL::Point_3<Mesh_K>>(p2, outFile);
				// Write attribute
				char buf[2] = {0, 0};
				outFile.write(buf, 2); // Set to zero
			}
			outFile.close();
		}

		void plotMeshMatlab2d(const TOMesh* const inMesh, std::ofstream& plotFile)
		{
			plotFile << "points = [";
			for(std::size_t k = 0; k < inMesh->getNumNodes(); ++k)
			{
				Mesh_K::Point_2 p0 = inMesh->getNode2D(k);
				plotFile << std::setprecision(15) << p0.x() << " " << p0.y() << '\n';
			}
			plotFile << "];" << '\n';
			plotFile << "tris = [";
			for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			{
				std::vector<std::size_t> elems = inMesh->getElementConnectivity(k);
				for(std::size_t kn = 0; kn < elems.size(); ++kn )
					plotFile << elems[kn] + 1 << " ";
				plotFile << '\n';
			}
			plotFile << "];" << '\n';
			plotFile << "mats = [";
			for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			{
				plotFile << inMesh->getOptVal(k) << '\n';
			}
			plotFile << "];" << '\n';
			plotFile << "figure" << '\n';
			plotFile << "for k = 1:length(tris)" << '\n';
			plotFile << "	patch(points(tris(k,:),1), points(tris(k,:),2),mats(k));" << '\n';
			plotFile << "end" << '\n';
			plotFile << "axis equal" << '\n';
		}
	
		void plotMeshMatlab3d(const TOMesh* const inMesh, std::ofstream& plotFile)
		{
			plotFile << "points = [";
			for(std::size_t k = 0; k < inMesh->getNumNodes(); ++k)
			{
				Mesh_K::Point_3 p0 = inMesh->getNode3D(k);
				plotFile << std::setprecision(15) << p0.x() << " " << p0.y() << " " << p0.z() << '\n';
			}
			plotFile << "];" << '\n';
			plotFile << "tets = [";
			for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			{
				std::vector<std::size_t> elems = inMesh->getElementConnectivity(k);
//				if(inMesh->getOptVal(k) != 0)
				{
					for(std::size_t kn = 0; kn < elems.size(); ++kn )
						plotFile << elems[kn] + 1 << " ";
					plotFile << '\n';
				}
			}
			plotFile << "];" << '\n';
			plotFile << "mats = [";
			for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			{
//				if(inMesh->getOptVal(k) != 0)
					plotFile << inMesh->getOptVal(k) << '\n';
			}
      plotFile << "];" << '\n';
			plotFile << "figure;" << '\n';
			plotFile << "caxis([0. 1.]);" << '\n';
			plotFile << "tetramesh(tets, points, mats);" << '\n';
		}

		unsigned getVTKCellType(unsigned dim, std::size_t nnodes)
		{
			if(dim == 2)
			{
				if(nnodes == 3)
					return 5; // Triangle
				else if(nnodes == 4)
					return 9; // Quadrilateral
			}
			if(nnodes == 4)
				return 10; // Tetrahedron
			else if(nnodes == 8)
				return 12; // Hexahedron
			return 0; // Not found
		}
		
		unsigned getGMSHCellType(unsigned dim, std::size_t nnodes)
		{
			if(dim == 2)
      {
        if(nnodes == 3)
          return 2; // Triangle
        else if(nnodes == 4)
          return 3; // Quadrilateral
      }
      if(nnodes == 4)
        return 4; // Tetrahedron
      else if(nnodes == 8)
        return 5; // Hexahedron
      return 0; // Not found
		}
	}

	void plotMeshMatlab(const TOMesh* const inMesh, const std::string& fileName)
	{
		std::ofstream plotFile(fileName.c_str());
		plotFile.precision(15);
		plotMeshMatlab(inMesh, plotFile);
		plotFile.close();
	}

	void plotMeshMatlab(const TOMesh* const inMesh, std::ofstream& plotFile)
	{
		if(inMesh->dimNum() == 2)
			plotMeshMatlab2d(inMesh, plotFile);
		else if(inMesh->dimNum() == 3)
			plotMeshMatlab3d(inMesh, plotFile);
	}

	void plotSegmentsMatlab(const std::vector<std::vector<Mesh_Segment_2> >& segVV, const std::string& fileName)
	{
		char plotColors[7] = {'r', 'b', 'g', 'm', 'c', 'y', 'k'};
		std::ofstream outFile(fileName.c_str());
		outFile << "figure;" << '\n';
		outFile << "hold on" << '\n';
		unsigned colorCount = 0, numColors = 7;
		for(Uint k1 = 0; k1 < segVV.size(); ++k1)
		{
			outFile << "segs" << k1 <<" = [";
			for(Uint k2 = 0; k2 < segVV[k1].size(); ++k2)
			{
				Point_2_base p1 = segVV[k1][k2].source(), p2 = segVV[k1][k2].target();
				outFile << std::setprecision(15) << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << '\n';
			}
			outFile << "];" << '\n';
			outFile << "for k = 1:size(segs" << k1 <<", 1)" << '\n';
			outFile << "  plot(segs" << k1 << "(k, [1 3]), segs" << k1 << "(k, [2 4]), '" << plotColors[colorCount] << "');" << '\n';
			outFile << "end" << '\n';
			++colorCount;
			if(colorCount >= numColors)
				colorCount = 0;
		}
		outFile << "hold off" << '\n';
		outFile << "axis equal" << '\n';
	}

	void plotSegmentsMatlab(const std::vector<Mesh_Segment_2>& segVec, const std::string& fileName)
	{
		std::ofstream outFile(fileName.c_str());
	  outFile << "segs = [";
	  for(Uint k = 0; k < segVec.size(); k++)
	  {
		Point_2_base p1 = segVec[k].source(), p2 = segVec[k].target();
		outFile << std::setprecision(15) << p1.x() << " " << p1.y() << " " << p2.x() << " " << p2.y() << '\n';
	  }
	  outFile << "];" << '\n';
	  outFile << "figure;" << '\n';
	  outFile << "hold on" << '\n';
	  outFile << "for k = 1:size(segs, 1)" << '\n';
	  outFile << "  plot(segs(k, [1 3]), segs(k, [2 4]))" << '\n';
	  outFile << "end" << '\n';
	  outFile << "hold off" << '\n';
	  outFile << "axis equal" << '\n';
	}

	void plotDataMatlab(const std::vector<double>& data, std::size_t nx, std::size_t ny, std::size_t nz, const std::string& fileName)
	{
		std::ofstream outFile(fileName);
		outFile.precision(15);
		outFile << "data = zeros(" << ny << ", " << nx << ", " << nz << ");" << '\n';
		for(std::size_t kz = 0; kz < nz; ++kz)
		{
			outFile << "data(:,:," << kz+1 << ") = [";
			for(std::size_t ky = 0; ky < ny; ++ky)
			{
				for(std::size_t kx = 0; kx < nx; ++kx)
				{
					std::size_t kp = kz + nz*ky + nz*ny*kx;
					outFile << data[kp] << " ";
				}
				outFile << '\n';
			}
			outFile << "];" << '\n';
		}
		outFile.close();
	}

	void plotDataMatlab(const std::vector<double>& data, std::size_t nx, std::size_t ny, const std::string& fileName)
	{
		std::ofstream outFile(fileName);
		outFile.precision(15);
		outFile << "data = [";
		for(std::size_t ky = 0; ky < ny; ++ky)
		{
			for(std::size_t kx = 0; kx < nx; ++kx)
			{
				std::size_t kp = ky + ny*kx;
				outFile << data[kp] << " ";
			}
			outFile << '\n';
		}
		outFile << "];" << '\n';
		outFile.close();
	}

	void plotSTL(const TOMesh* const inMesh, const std::string& fileName, bool useASCII)
	{
		plotSTL(*inMesh, fileName, useASCII);
	}

	void plotSTL(const TOMesh& inMesh, const std::string& fileName, bool useASCII)
	{
		if(useASCII)
			plotSTLASCII(inMesh, fileName);
		else
			plotSTLBinary(inMesh, fileName);
	}

	void plotSTL(const C2t3& surfMesh, const std::string& fileName, bool useASCII)
	{
		if(useASCII)
			plotSTLASCII(surfMesh, fileName);
		else
			plotSTLBinary(surfMesh, fileName);
	}
	
	void plotMeshVTK(const TOMesh* const inMesh, unsigned kt, const std::string& outFileName)
	{
		std::stringstream ss;
		ss << outFileName << "." << kt << ".vtk";
		plotMeshVTK(inMesh, ss.str());
	}

	void plotMeshVTK(const TOMesh* const inMesh, const std::string& outFileName)
	{
		std::ofstream outFile(outFileName.c_str());
		plotMeshVTK(inMesh, outFile);
	}

	void plotMeshVTK(const TOMesh* const inMesh, std::ofstream& outFile)
	{
		outFile.precision(8);
		outFile << std::scientific;
		outFile << "# vtk DataFile Version 3.1" << '\n';
		outFile << "Output for topologies\nASCII" << '\n';
		outFile << "DATASET UNSTRUCTURED_GRID" << '\n';
		outFile << "POINTS " << inMesh->getNumNodes() << " FLOAT" << '\n';
		// Output points
		for(std::size_t k = 0; k < inMesh->getNumNodes(); ++k)
		{
			if(inMesh->dimNum() == 2)
				outFile << inMesh->getNode2D(k) << " 0." << '\n';
			else
				outFile << inMesh->getNode3D(k) << '\n';
		}
		// Output elements
		// First count number of cell nodes
		std::size_t nnodes = 0;
		for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			nnodes += inMesh->getElementConnectivity(k).size();
		std::size_t ndata = nnodes + inMesh->getNumElements(); // Number of numbers
		// Ouput cell information
		outFile << "CELLS " << inMesh->getNumElements() << " " << ndata << '\n';
		for(std::size_t ke = 0; ke < inMesh->getNumElements(); ++ke)
		{
			const std::vector<std::size_t>& curCell = inMesh->getElementConnectivity(ke);
			outFile << curCell.size();
			for(std::size_t kc = 0; kc < curCell.size(); ++kc)
				outFile << " " << curCell[kc];
			outFile << '\n';
		}
		// Output cell types, currently only tris and tets are supported
		outFile << "CELL_TYPES " << inMesh->getNumElements() << '\n';
		for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			outFile << getVTKCellType(inMesh->dimNum(), inMesh->getElementConnectivity(k).size()) << '\n';
		// Output mesh data
		outFile << "CELL_DATA " << inMesh->getNumElements() << '\n';
		// Output headers
		outFile << "SCALARS density FLOAT 1" << '\n';
		outFile << "LOOKUP_TABLE default" << '\n';
		// Output data
		for(std::size_t ke = 0; ke < inMesh->getNumElements(); ++ke)
			outFile << inMesh->getOptVal(ke) << '\n';
	}

	void plotMeshGMSH(const TOMesh* const inMesh, const std::string& outFileName, unsigned kt, bool writeMesh)
	{
		std::ios_base::openmode mode = writeMesh ? std::ios_base::out : (std::ios_base::app | std::ios_base::out);
		std::ofstream outFile(outFileName.c_str(), mode);
		plotMeshGMSH(inMesh, outFile, kt, writeMesh);
	}

	void plotMeshGMSH(const TOMesh* const inMesh, std::ofstream& outFile, unsigned kt, bool writeMesh)
	{
		outFile.precision(8);
		outFile << std::scientific;
		if(writeMesh) // Write mesh
		{
			outFile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat" << '\n'; // Header info
			// Write nodes
			outFile << "$Nodes\n" << inMesh->getNumNodes() << '\n';
			for(std::size_t k = 0; k < inMesh->getNumNodes(); ++k)
			{
				outFile << k+1 << " ";
				if(inMesh->dimNum() == 2)
					outFile << inMesh->getNode2D(k) << " 0." << '\n';
				else
					outFile << inMesh->getNode3D(k) << '\n';
			}
			outFile << "$EndNodes" << '\n';
			// Write elements
			outFile << "$Elements\n" << inMesh->getNumElements() << '\n';
			for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			{
				const std::vector<std::size_t>& curElem = inMesh->getElementConnectivity(k);
				outFile << k+1 << " " << getGMSHCellType(inMesh->dimNum(), curElem.size()) << " ";
				outFile << "2 " << inMesh->getMatID(k) << " " << inMesh->getMatID(k);
				for(std::size_t kc = 0; kc < curElem.size(); ++kc)
					outFile << " " << curElem[kc]+1;
				outFile << '\n';
			}
			outFile << "$EndElements"<< '\n';
		}
		// Write data
		outFile << "$ElementData\n1\nDensity" << '\n';
		outFile << "1\n" << kt << '\n';
		outFile << "3\n" << kt << "\n1" << '\n';
		outFile << inMesh->getNumElements() << '\n';
		for(std::size_t k = 0; k < inMesh->getNumElements(); ++k)
			outFile << k+1 << " " << inMesh->getOptVal(k) << '\n';
		outFile << "$EndElementData" << '\n';
	}

	void outputTORData(const TopOptRep* outTOR, const std::string& fileName)
	{
		std::ofstream outFile(fileName);
    outFile.precision(16);
		std::vector<std::vector<int> > discreteVars;
		std::vector<std::vector<double> > realVars;
		outTOR->getMPIRep(discreteVars, realVars);
		outFile << discreteVars.size() << '\n';
		for(std::size_t k1 = 0; k1 < discreteVars.size(); ++k1)
		{
			outFile << discreteVars[k1].size() << '\n';
			for(std::size_t k2 = 0; k2 < discreteVars[k1].size(); ++k2)
				outFile << discreteVars[k1][k2] << " ";
			outFile << '\n';
		}
		outFile << realVars.size() << '\n';
    for(std::size_t k1 = 0; k1 < realVars.size(); ++k1)
    {
      outFile << realVars[k1].size() << '\n';
      for(std::size_t k2 = 0; k2 < realVars[k1].size(); ++k2)
        outFile << realVars[k1][k2] << " ";
      outFile << '\n';
    }
	}
}
}

