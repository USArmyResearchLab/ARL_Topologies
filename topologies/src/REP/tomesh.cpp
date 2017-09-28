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

#include "tomesh.h"
#include "helper.h"
#include <map>

namespace Topologies{
// TOMesh base class
TOMesh::TOMesh(unsigned inDimNum, std::vector<std::vector<std::size_t>> inCVec) : 
	m_dimNum(inDimNum), 
	connectivityVec(std::move(inCVec)), 
	optVals(std::vector<double>(connectivityVec.size(), 0.)), 
	matIDVec(std::vector<int>(connectivityVec.size(), 0)) 
{
}

TOMesh::TOMesh(unsigned inDimNum, std::vector<std::vector<std::size_t>> inCVec, std::vector<double> inOVec) : 
	m_dimNum(inDimNum), 
	connectivityVec(std::move(inCVec)), 
	optVals(std::move(inOVec)), 
	matIDVec(std::vector<int>(connectivityVec.size(), 0)) 
{
}

TOMesh::TOMesh(unsigned inDimNum, std::vector<std::vector<std::size_t>> inCVec, std::vector<int> inMIDVec) :
	m_dimNum(inDimNum), 
	connectivityVec(std::move(inCVec)), 
	optVals(std::vector<double>(connectivityVec.size(), 0.)), 
	matIDVec(std::move(inMIDVec)) 
{
}

TOMesh::TOMesh(unsigned inDimNum, std::vector<std::vector<std::size_t>> inCVec, std::vector<double> inOVec, 
					std::vector<int> inMIDVec) :
	m_dimNum(inDimNum), 
	connectivityVec(std::move(inCVec)), 
	optVals(std::move(inOVec)), 
	matIDVec(std::move(inMIDVec)) 
{
}

void TOMesh::setOptVals(const std::vector<double>& inVec)
{
	if(inVec.size() == optVals.size())
		optVals = inVec;
}

void TOMesh::TOMesh::setOptVal(std::size_t ke, double val)
{
	assert(ke < optVals.size());
	optVals[ke] = val;
}

void TOMesh::setMatIDs(const std::vector<int>& inVec)
{
	if(inVec.size() == matIDVec.size())
		matIDVec = inVec;
}

void TOMesh::setMatIDs(std::size_t ke, int val)
{
	assert(ke < matIDVec.size());
	matIDVec[ke] = val;
}

void TOMesh::swap(TOMesh& arg)
{
	connectivityVec.swap(arg.connectivityVec);
	optVals.swap(arg.optVals);
	matIDVec.swap(arg.matIDVec);
	std::swap(m_dimNum, arg.m_dimNum);
}

// TOMesh2D
TOMesh2D::TOMesh2D(const CDT_2& inMesh) : TOMesh(2) 
{
	setupFromCDT2(inMesh);
	matIDVec = std::vector<int>(optVals.size(), 0);
}

TOMesh2D::TOMesh2D(std::vector<Point_2_base> inPtVec, std::vector<std::vector<std::size_t>> inConnVec) : 
	TOMesh(2, std::move(inConnVec)), 
	ptVec(std::move(inPtVec)) 
{
}

TOMesh2D::TOMesh2D(std::vector<Point_2_base> inPtVec, std::vector<std::vector<std::size_t>> inConnVec, 
			std::vector<double> inOptVec) :
	TOMesh(2, std::move(inConnVec), std::move(inOptVec)), 
	ptVec(std::move(inPtVec)) 
{
}

TOMesh2D::TOMesh2D(std::vector<Point_2_base> inPtVec, std::vector<std::vector<std::size_t>> inConnVec,
		std::vector<int> inMIDVec) :
	TOMesh(2, std::move(inConnVec), std::move(inMIDVec)), 
	ptVec(std::move(inPtVec)) 
{
}

TOMesh2D::TOMesh2D(std::vector<Point_2_base> inPtVec, std::vector<std::vector<std::size_t>> inConnVec,
		std::vector<double> inOptVec, std::vector<int> inMIDVec) :
	TOMesh(2, std::move(inConnVec), std::move(inOptVec), std::move(inMIDVec)), 
	ptVec(std::move(inPtVec)) 
{
}

void TOMesh2D::swap(TOMesh2D& arg)
{
	TOMesh::swap(arg);
	ptVec.swap(arg.ptVec);
}

void TOMesh2D::setupFromCDT2(const CDT_2& inMesh)
{
	// Store points
	ptVec.reserve(inMesh.number_of_vertices());
	for(CDT_2::Point_iterator pit = inMesh.points_begin(); pit != inMesh.points_end(); ++pit)
		ptVec.push_back(*pit);
	// Store elements and opt. vals.
	connectivityVec.reserve(inMesh.number_of_faces());
	optVals.reserve(inMesh.number_of_faces());
	for(CDT_2::Finite_faces_iterator fit = inMesh.finite_faces_begin(); fit != inMesh.finite_faces_end(); ++fit)
	{
		if(fit->is_in_domain())
		{
			std::vector<std::size_t> nodeIds(3);
			for(unsigned k = 0; k < 3; ++k)
				nodeIds[k] = fit->vertex(k)->info();
			connectivityVec.push_back(nodeIds);
			optVals.push_back(fit->info().optVal);
		}
	}
}

// TOMesh3D
TOMesh3D::TOMesh3D(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems) :
	TOMesh(3, std::move(inElems)), 
	ptVec(std::move(inPts)) 
{
}
TOMesh3D::TOMesh3D(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems, 
	std::vector<double> inOptVals) :
	TOMesh(3, std::move(inElems), std::move(inOptVals)), 
	ptVec(std::move(inPts))	
{
}
TOMesh3D::TOMesh3D(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems,
		std::vector<int> inMIDVec) :
	TOMesh(3, std::move(inElems), std::move(inMIDVec)), 
	ptVec(std::move(inPts)) 
{
}

TOMesh3D::TOMesh3D(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems,
		std::vector<double> inOptVals, std::vector<int> inMIDVec) :
	TOMesh(3, std::move(inElems), std::move(inOptVals), std::move(inMIDVec)), 
	ptVec(std::move(inPts))
{
}

TOMesh3D::TOMesh3D(const C3t3& inMesh) : TOMesh(3)
{
	setupFromC3t3(inMesh);
	matIDVec = std::vector<int>(connectivityVec.size(), 0);
	optVals = std::vector<double>(connectivityVec.size(), 0.);
}

void TOMesh3D::swap(TOMesh3D& arg)
{
	TOMesh::swap(arg);
	ptVec.swap(arg.ptVec);
	edgeVec.swap(arg.edgeVec);
	faceVec.swap(arg.faceVec);
}

void TOMesh3D::setupFromC3t3(const C3t3& inMesh)
{
	const Tr& tr = inMesh.triangulation();
	ptVec.reserve(tr.number_of_vertices());
	std::map<Tr::Vertex_handle, std::size_t> handleIdMap;
	std::size_t count = 0;
	// Copy points and get node ids
	for(auto vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit)
	{
		ptVec.push_back(vit->point());
		handleIdMap[vit] = count++;
	}
	// Set up tetrahedra
	connectivityVec.reserve(inMesh.number_of_cells_in_complex());
	for(auto cit = inMesh.cells_in_complex_begin(); cit != inMesh.cells_in_complex_end(); ++cit)
		connectivityVec.push_back({handleIdMap[cit->vertex(0)], handleIdMap[cit->vertex(1)], 
															handleIdMap[cit->vertex(2)], handleIdMap[cit->vertex(3)]});
	
}

// TOMesh3DSurface
TOMesh3DSurface::TOMesh3DSurface(const C2t3& inMesh) : TOMesh(3)
{
	setupFromC2t3(inMesh);
}

TOMesh3DSurface::TOMesh3DSurface(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems):
	TOMesh(3, std::move(inElems)),
	ptVec(std::move(inPts))
{
}

void TOMesh3DSurface::swap(TOMesh3DSurface& arg)
{
	TOMesh::swap(arg);
	ptVec.swap(arg.ptVec);
}

void TOMesh3DSurface::setupFromC2t3(const C2t3& inMesh)
{
	const TOM_Tr& tr = inMesh.triangulation();
	std::map<const TOM_Vertex_handle, std::size_t> V;
	ptVec.reserve(tr.number_of_vertices());
	std::size_t i = 0;
	for(TOM_Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit)
	{
		ptVec.push_back(vit->point());
		V[vit] = i++;
	}
	for(TOM_Tr::Finite_facets_iterator fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
	{
		const typename TOM_Tr::Cell_handle cell = fit->first;
		const int& index = fit->second;
		std::vector<std::size_t> indexVec(3);
		if(cell->is_facet_on_surface(index)==true)
		{
			for(unsigned k = 0; k < 3; ++k)
			{
				const TOM_Vertex_handle& tvh = cell->vertex(tr.vertex_triple_index(index, k));
				std::map<const TOM_Vertex_handle, std::size_t>::const_iterator mit = V.find(tvh);
				if(mit != V.end())
					indexVec[k] = mit->second;
			}
			connectivityVec.push_back(indexVec);
		}
		//			connectivityVec.push_back(indexVec);
	}
	optVals = std::vector<double>(connectivityVec.size(), 0.);
	matIDVec = std::vector<int>(connectivityVec.size(), 0);
}
}// namespace

