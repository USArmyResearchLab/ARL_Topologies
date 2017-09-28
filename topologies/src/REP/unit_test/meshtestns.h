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

#ifndef MESHTESTNS_H
#define MESHTESTNS_H

#include "tomesh.h"
#include "tomeshprocessing.h"
#include "inputloaderrep.h"
#include "cgal_types.h"

namespace MeshTestNS
{
	using namespace Topologies; 

	InputLoader::RepNodeInfo loadRNI(const std::string& fileName)
	{
		// Loads an XML document and parses the first <representation> node seen
		pugi::xml_document xmldoc;
		xmldoc.load_file(fileName.c_str());
		assert(xmldoc);
		pugi::xml_node rootNode = xmldoc.child("representation");
		InputLoader::RepNodeInfo testRNI;
		testRNI.parse(rootNode, fileName);
		return testRNI;
	}

	std::pair<double,double> getXRange(const CDT_2& mesh)
	{
		double maxx, minx;
		bool first = true;
		for(auto it = mesh.points_begin(); it != mesh.points_end(); ++it)
		{
			if(first || (it->x() > maxx))
				maxx = it->x();
			if(first || (it->x() < minx))
				minx = it->x();
			first = false;
		}
		return std::make_pair(minx, maxx);
	}

	std::pair<double,double> getYRange(const CDT_2& mesh)
	{
		double maxy, miny;
		bool first = true;
		for(auto it = mesh.points_begin(); it != mesh.points_end(); ++it)
		{
			if(first || (it->y() > maxy))
				maxy = it->y();
			if(first || (it->y() < miny))
				miny = it->y();
			first = false;
		}
		return std::make_pair(miny, maxy);
	}

	std::vector<double> getBoundBox(const C3t3& mesh)
	{
		const Tr& tr = mesh.triangulation();
		double minx, maxx, miny, maxy, minz, maxz;
		bool first = true;
		for(auto vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit)
		{
			double x = vit->point().x(), y = vit->point().y(), z = vit->point().z();
			if(first || (x < minx))
				minx = x;
			if(first || (x > maxx))
				maxx = x;
			if(first || (y < miny))
				miny = y;
			if(first || (y > maxy))
				maxy = y;
			if(first || (z < minz))
				minz = z;
			if(first || (z > maxz))
				maxz = z;
			first = false;
		}
		return {minx, maxx, miny, maxy, minz, maxz};
	}

	double getMeshArea(const CDT_2& mesh)
	{
		double area = 0.;
		for(auto it = mesh.finite_faces_begin(); it != mesh.finite_faces_end(); ++it)
			if(it->is_in_domain())
				area += CGAL::area(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point());
		return area;
	}

	double getMeshArea(TOMesh& inMesh)
	{
		double area = 0.;
		for(std::size_t k = 0; k < inMesh.getNumElements(); ++k)
			area += TOMeshProcessing::computeElementVolume(k, &inMesh);
		return area;
	}

	double getMeshVolume(const C3t3& mesh)
	{
		double vol = 0.;
		for(auto cit = mesh.cells_in_complex_begin(); cit != mesh.cells_in_complex_end(); ++cit)
			vol += CGAL::volume<CGAL::Epick>(cit->vertex(0)->point(), cit->vertex(1)->point(), cit->vertex(2)->point(), cit->vertex(3)->point());
		return vol;
	}

	double computeMaterialPropertyTotal(const CDT_2& mesh, unsigned k, unsigned numProps)
	{
		double propTotal = 0.;
		for(auto it = mesh.finite_faces_begin(); it != mesh.finite_faces_end(); ++it)
		{
			const GenericMaterial& curGM = it->info().triMat;
			assert(curGM.getNumParameters() == numProps);
			double area = CGAL::area(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point());
			propTotal += curGM.getParameter(k)*area;
		}
		return propTotal;
	}

	double computeOptValTotal(const CDT_2& mesh)
	{
		double propTotal = 0.;
		for(auto it = mesh.finite_faces_begin(); it != mesh.finite_faces_end(); ++it)
		{
			double area = CGAL::area(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point());
			propTotal += it->info().optVal*area;
		}
		return propTotal;
	}

	double computeSurfaceArea(const TOMesh3DSurface& tom3s)
	{
		double sa = 0.;
		for(std::size_t k = 0; k < tom3s.getNumElements(); ++k)
		{
			const std::vector<std::size_t>& connVec = tom3s.getElementConnectivity(k);
			assert(connVec.size() == 3);
			Mesh_K::Triangle_3 t3(tom3s.getNode3D(connVec[0]), tom3s.getNode3D(connVec[1]), tom3s.getNode3D(connVec[2]));
			sa += sqrt(t3.squared_area());
		}
		return sa;
	}

	void printTOMesh(const TOMesh& mesh, bool is3d = false)
	{
		std::cout << "pts = [";
		for(std::size_t k = 0; k < mesh.getNumNodes(); ++k)
		{
			if(is3d)
				std::cout << mesh.getNode3D(k) << std::endl;
			else
				std::cout << mesh.getNode2D(k) << std::endl;
		}
		std::cout << "];" << std::endl;
		std::cout << "elems = [";
		for(std::size_t ke = 0; ke < mesh.getNumElements(); ++ke)
		{
			std::vector<std::size_t> curElem = mesh.getElementConnectivity(ke);
			for(std::size_t kv = 0; kv < curElem.size(); ++kv)
				std::cout << curElem[kv] << " ";
			std::cout << mesh.getOptVal(ke) << std::endl;
		}
		std::cout << "];" << std::endl;
	}
}

#endif

