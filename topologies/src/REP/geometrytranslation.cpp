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

#include "geometrytranslation.h"
#include <list>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/intersections.h>

using std::vector;
using std::cout;
using std::endl;

namespace Topologies{
namespace GeometryTranslation
{

namespace // Anonymous namespace for helper functions
{
	std::vector<Mesh_Segment_2> getSegRep(unsigned nx, unsigned ny, double width, double height)
	{
		// Constraints are all box sides, just need to make sure we don't insert duplicates
		// Move along horizontal segments first
		std::vector<Mesh_Segment_2> segVec;
		segVec.reserve(nx*(ny + 1) + (nx + 1)*ny);
		double dx = width/(double)nx, dy = height/(double)ny;
		for(unsigned ky = 0; ky < (ny + 1); ++ky)
		{
			double cury = dy*(double)ky;
			for(unsigned kx = 0; kx < nx; ++kx)
			{
				double curx = dx*(double)kx;
				Point_2_base p1(curx, cury), p2(curx + dx, cury);
				segVec.push_back(Mesh_Segment_2(p1, p2));
			}
		}
		// Vertical segments
		for(unsigned kx = 0; kx < (nx + 1); ++kx)
		{
			double curx = dx*(double)kx;
			for(unsigned ky = 0; ky < ny; ++ky)
			{
				double cury = dy*(double)ky;
				Point_2_base p1(curx, cury), p2(curx, cury + dy);
				segVec.push_back(Mesh_Segment_2(p1, p2));
			}
		}
		return segVec;
	}

	vector<Point_2> orderedMeshSegToPointVec(const vector<Mesh_Segment_2>& segVec)
	{
		vector<Point_2> ptVec;
		for(unsigned k = 0; k < segVec.size(); ++k)
		{
			Point_2 curpt(segVec[k].source().x(), segVec[k].source().y());
			ptVec.push_back(curpt);
		}
		return ptVec;
	}
	
	bool pointInVector(const Point_2_base& chkPt, const vector<Point_2_base>& ptVec, unsigned& vid)
	{
		for(unsigned k = 0; k < ptVec.size(); ++k)
		{
			if(TOL_EQ(chkPt.x(), ptVec[k].x(), 1e-8) && TOL_EQ(chkPt.y(), ptVec[k].y(), 1e-8))
			{
				vid = k;
				return true;
			}
		}
		return false;
	} 
	
	unsigned addToPointVec(const Point_2_base& newPt, vector<Point_2_base>& ptVec, vector<unsigned>& countVec)
	{
		unsigned foundID = 0;
    if(!pointInVector(newPt, ptVec, foundID))
    {
      ptVec.push_back(newPt);
      countVec.push_back(1);
			foundID = countVec.size() - 1;
    }
    else
      countVec[foundID]++;
		return foundID;
	}
	
	CDT_2::Vertex_handle insertVertexInCDT2(CDT_2& inCDT, CDT_2::Point& p, const double bsnap2)
	{
		CDT_2::Vertex_handle vcur;
		CDT_2::Finite_vertices_iterator fvit;
		bool found = false;
		for(fvit = inCDT.finite_vertices_begin(); fvit != inCDT.finite_vertices_end() && !found; ++fvit)
		{
			if((fvit->point() - p).squared_length() < bsnap2) // Don't insert repeated vertices (crashes mesher)
			{
				vcur = fvit;
				found = true;
			}
		}
		if(!found)
			vcur = inCDT.insert(p);
		return vcur;
	}
	
	void mark_domains(CDT_2& ct, CDT_2::Face_handle start, int index, std::list<CDT_2::Edge>& border)
	{
		if(start->info().nesting_level != -1)
			return;
		std::list<CDT_2::Face_handle> queue;
		queue.push_back(start);
		while(! queue.empty())
		{
			CDT_2::Face_handle fh = queue.front();
			queue.pop_front();
			if(fh->info().nesting_level == -1)
			{
				fh->info().nesting_level = index;
				for(int i = 0; i < 3; i++)
				{
					CDT_2::Edge e(fh,i);
					CDT_2::Face_handle n = fh->neighbor(i);
					if(n->info().nesting_level == -1)
					{
						if(ct.is_constrained(e))
							border.push_back(e);
						else
							queue.push_back(n);
					}
				}
			}
		}
	}

	void mark_domains(CDT_2& cdt)
	{
		std::list<FaceInfo2> meshFI2;
		// Initialize nesting level
		for(CDT_2::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it)
			it->info().nesting_level = -1;
		int index = 0;
		std::list<CDT_2::Edge> border;
		// Mark infinite faces with 0, not to be used in final mesh
		mark_domains(cdt, cdt.infinite_face(), index++, border);
		// Loop through and mark remaining faces
		while(! border.empty())
		{
			CDT_2::Edge e = border.front();
			border.pop_front();
			CDT_2::Face_handle n = e.first->neighbor(e.second);
			if(n->info().nesting_level == -1)
				mark_domains(cdt, n, e.first->info().nesting_level+1, border);
		}
	}

	void setTriangleInfo(CDT_2& cdt, const GenericMaterial& inGM1, const GenericMaterial inGM2) // 2 materials, holes get 2nd
	{
		for(CDT_2::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it)
		{
			it->info().optVal = 1.;
			if(it->info().nesting_level == 1)
				it->info().triMat = inGM1;
			else 
				it->info().triMat = inGM2;
		}
	}
	void setTriangleInfo(CDT_2& cdt, const GenericMaterial& inGM1) // Sets all triangles to 1 material
	{
		for(CDT_2::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it)
		{
			it->info().optVal = 1.;
			it->info().triMat = inGM1;
		}
	}

	CDT_2 mesh2D(const std::vector<Mesh_Segment_2>& inSegVec, std::list<CDT_2::Point> seedPts, 
							 const MesherData& meshParams)
	{
		bool debugMesh = false;
		double bsnap = meshParams.triMeshEdgeSize*0.05; // Minimum constraint size to prevent small triangles
		double bsnap2 = bsnap*bsnap;
		if(debugMesh)
		{
			cout << "insert points and constraints into the CDT" << endl;
			cout << "Number of shapes: " << inSegVec.size() << endl;
			cout << "figure; hold on;" << endl;
		}
		bool usered = false;
		CDT_2 outMesh;
		for(unsigned k = 0; k < inSegVec.size(); ++k)
		{
			Mesh_Segment_2 curseg = inSegVec[k];
			CDT_2::Point p1(curseg.source());
			CDT_2::Vertex_handle v1 = insertVertexInCDT2(outMesh, p1, bsnap2);
			CDT_2::Point p2(curseg.target());
			CDT_2::Vertex_handle v2 = insertVertexInCDT2(outMesh, p2, bsnap2);
			p1 = v1->point();
			p2 = v2->point();
			Real segLen2 = (p1 - p2).squared_length();
			if(segLen2 >= bsnap2)
			{
				outMesh.insert_constraint(v1, v2);
				if(debugMesh)
				{
					cout << "p1 =[" << p1 << "];" << endl;
					cout << "p2 = [" << p2 << "];" << endl;
					if(usered)
						cout << "plot([p1(1) p2(1)],[p1(2) p2(2)],'r.-')" << endl;
					else
						cout << "plot([p1(1) p2(1)],[p1(2) p2(2)],'b.-')" << endl;
					usered = !usered;
					cout << "pause" << endl;
				}
			}
		}
		if(debugMesh) cout << "Generating mesh" << endl;
		Mesher_2 mesher(outMesh);
		Real b = 2.*pow(sin(Deg2Rad*meshParams.triMeshEdgeAngle), 2.);
		if(debugMesh)
			cout << "edge size: " << meshParams.triMeshEdgeSize << endl;
		mesher.set_criteria(Criteria(b, meshParams.triMeshEdgeSize));
		if(!seedPts.empty())
			mesher.set_seeds(seedPts.begin(), seedPts.end());
		if(debugMesh)
		{
			cout << "Mesh before refinement: " << endl;
			cout << outMesh << endl;
			cout << "refining" << endl;
		}
		mesher.refine_mesh();

		// Default triangle info:
		for(CDT_2::Finite_faces_iterator it = outMesh.finite_faces_begin(); it != outMesh.finite_faces_end(); ++it)
			it->info().optVal = 1.;
		// Set vertex ids
		unsigned id = 0;
		for(CDT_2::Finite_vertices_iterator it = outMesh.finite_vertices_begin(); it != outMesh.finite_vertices_end(); ++it)
		{
			it->info() = id;
			id++;
		}
		return outMesh;
	}
	
	void computeSeedPoints(const std::vector<std::vector<Mesh_Segment_2> >& inSegVec, std::list<CDT_2::Point>& outPts)
	{
		if(inSegVec.size() == 1)
			return;
		// This assumes that the set of mesh segments in inSegVec consist of one outer shell and the remaining vectors 
		// are holes.  This will determine the outer shell and only compute seeds for the holes.
		std::vector<double> areas;
		for(unsigned k = 0; k < inSegVec.size(); ++k)
		{
			// Partition the polygon and use the centroid of a partition
			CGAL::Polygon_2<Mesh_K> polygon;
			meshSegs2Polygon(inSegVec[k], polygon);
			areas.push_back(fabs(polygon.area()));
			std::list<CGAL::Partition_traits_2<Mesh_K>::Polygon_2> partitionPolys;
			CGAL::approx_convex_partition_2(polygon.vertices_begin(), polygon.vertices_end(),
																			std::back_inserter(partitionPolys));
			if(!partitionPolys.empty())
			{
				CGAL::Partition_traits_2<Mesh_K>::Polygon_2 curPoly = *partitionPolys.begin();
				CGAL::Partition_traits_2<Mesh_K>::Polygon_2::Vertex_iterator iv;
				Point_2_base partitionCentroid(0., 0.);
				double np = 0.;
				for(iv = curPoly.vertices_begin(); iv != curPoly.vertices_end(); ++iv)
				{
					Point_2_base curpt = *iv;
					partitionCentroid = partitionCentroid + (curpt - CGAL::ORIGIN);
					np += 1.;
				}
				partitionCentroid = CGAL::ORIGIN + (partitionCentroid - CGAL::ORIGIN)/np;
				outPts.push_back(partitionCentroid);
			}
		}
		// Delete the seed point associated with the outer shell
		if(!areas.empty())
		{
			std::vector<double>::const_iterator maxit = std::max_element(areas.begin(), areas.end());
			outPts.erase(std::next(outPts.begin(), maxit - areas.begin()));
		}
	}
	
	void pushBackPointUnique(const Point_2_base& p, std::vector<Point_2_base>& pointVec, double tol2)
	{
		bool found = false;
		for(unsigned k = 0; k < pointVec.size() && !found; ++k)
		{
			if(CGAL::squared_distance(p, pointVec[k]) < tol2)
				found = true;
		}
		if(!found)
			pointVec.push_back(p);
	}
	
	void polygon2MeshSegs(const CGAL::Partition_traits_2<Mesh_K>::Polygon_2& inPoly, std::vector<Mesh_Segment_2>& outSegVec)
	{
		outSegVec.resize(inPoly.size());
		CGAL::Partition_traits_2<Mesh_K>::Polygon_2::Vertex_iterator iv = inPoly.vertices_begin();
		for(unsigned k = 0; k < inPoly.size(); ++k)
		{
			Point_2_base p1 = *iv;
			++iv;
			if(iv == inPoly.vertices_end())
				iv = inPoly.vertices_begin();
			outSegVec[k] = Mesh_Segment_2(p1, *iv);
		}
	}
} // end anonymous namespace

void partitionMeshSegPolys(const std::vector<Mesh_Segment_2>& inSegVec, std::vector<std::vector<Mesh_Segment_2> >& outSegVec)
{
	typedef std::list<CGAL::Partition_traits_2<Mesh_K>::Polygon_2> Poly_list;
	std::vector<std::vector<Mesh_Segment_2> > orderedVecVec;
	orderMeshSegments(inSegVec, orderedVecVec);
	for(unsigned k = 0; k < orderedVecVec.size(); ++k)
	{
		CGAL::Polygon_2<Mesh_K> polygon;
		meshSegs2Polygon(orderedVecVec[k], polygon);
		Poly_list partitionPolys;
		CGAL::approx_convex_partition_2(polygon.vertices_begin(), polygon.vertices_end(),
										std::back_inserter(partitionPolys));
		for(Poly_list::const_iterator plit = partitionPolys.begin(); plit != partitionPolys.end(); ++plit)
		{
			std::vector<Mesh_Segment_2> curSegVec;
			polygon2MeshSegs(*plit, curSegVec);
			outSegVec.push_back(curSegVec);
		}
	}
}

void meshSegs2Polygon(const std::vector<Mesh_Segment_2>& inSegVec, CGAL::Polygon_2<Mesh_K>& outPoly)
{
	outPoly.clear();
	for(unsigned k = 0; k < inSegVec.size(); ++k)
		outPoly.push_back(inSegVec[k].source());
}

void meshSegs2PolygonWithHoles(const std::vector<Mesh_Segment_2>& inSegVec, CGAL::Polygon_2<Mesh_K>& outPoly, std::vector<CGAL::Polygon_2<Mesh_K> >& holeVec)
{
	holeVec.clear();
	std::vector<std::vector<Mesh_Segment_2> > orderedBoundVec;
	orderMeshSegments(inSegVec, orderedBoundVec);
	std::vector<CGAL::Polygon_2<Mesh_K> > boundPolyVec(orderedBoundVec.size());
	for(unsigned k = 0; k < orderedBoundVec.size(); ++k)
		meshSegs2Polygon(orderedBoundVec[k], boundPolyVec[k]);
	// Get enclosing polygon (max area)
	if(!boundPolyVec.empty())
	{
		double maxArea = fabs(boundPolyVec[0].area());
		unsigned maxk = 0;
		for(unsigned k = 1; k < boundPolyVec.size(); ++k)
		{
			double curArea = fabs(boundPolyVec[k].area());
			if(curArea > maxArea)
			{
				maxk = k;
				maxArea = curArea;
			}
		}
		outPoly = boundPolyVec[maxk];
		if(outPoly.is_clockwise_oriented())
			outPoly.reverse_orientation();
		// Set up any holes
		for(unsigned k = 0; k < boundPolyVec.size(); ++k)
		{
			if(k != maxk)
			{
				if(boundPolyVec[k].is_counterclockwise_oriented())
					boundPolyVec[k].reverse_orientation();
				holeVec.push_back(boundPolyVec[k]);
			}
		}
	}
}

void polygon2MeshSegs(const CGAL::Polygon_2<Mesh_K>& inPoly, std::vector<Mesh_Segment_2>& outSegVec)
{
	unsigned nsegs = inPoly.size();
	outSegVec.resize(nsegs);
	for(unsigned k = 0; k < nsegs; ++k)
	{
		unsigned kp1 = (k + 1) % nsegs;
		Mesh_Segment_2 tmpSeg(inPoly[k], inPoly[kp1]);
		outSegVec[k] = tmpSeg;
	}
}

void projectPointToPolygon(Point_2_base& chckPt, const CGAL::Polygon_2<Mesh_K>& poly, double tol)
{
	typedef CGAL::Polygon_2<Mesh_K>::Edge_const_iterator Edge_iterator;
	Edge_iterator closestEIT = poly.edges_begin();
	double closestDist = CGAL::squared_distance(chckPt, *poly.edges_begin());
	for(Edge_iterator eit = poly.edges_begin(); eit != poly.edges_end(); ++eit)
	{
		double curDist = CGAL::squared_distance(chckPt, *eit);
		if(curDist < closestDist)
		{
			closestEIT = eit;
			closestDist = curDist;
		}
	}
	chckPt = closestEIT->supporting_line().projection(chckPt);
	double chckDist2 = CGAL::squared_distance(chckPt, *closestEIT);
	if(chckDist2 > tol*tol)
	{
		// Still outside polygon so set to closest vertex of polygon
		Point_2_base p1 = closestEIT->source(), p2 = closestEIT->target();
		if(CGAL::has_smaller_distance_to_point(chckPt, p1, p2))
			chckPt = p1;
		else
			chckPt = p2;
	}
}

void polygon2ComponentVec(const CGAL::Polygon_2<Mesh_K>& inPoly, std::vector<double>& outComponentVec)
{
	outComponentVec.resize(2*inPoly.size());
  for(unsigned kp = 0; kp < inPoly.size(); ++kp)
  {
    outComponentVec[2*kp] = inPoly[kp].x();
    outComponentVec[2*kp + 1] = inPoly[kp].y();
  }
}

void componentVec2Polygon(const std::vector<double>& inComponentVec, CGAL::Polygon_2<Mesh_K>& outPoly)
{
	outPoly.clear();
	if((inComponentVec.size() % 2) == 0)
	{
		for(unsigned k = 0; k < inComponentVec.size(); k+=2)
			outPoly.push_back(Point_2_base(inComponentVec[k], inComponentVec[k+1]));
	}
}

Real computeSignedArea(const vector<Mesh_Segment_2>& inSegVec)
{
	Real a = 0.;
	vector<Mesh_Segment_2>::const_iterator cvit;
	for(cvit = inSegVec.begin(); cvit != inSegVec.end(); ++cvit)
	{
		Point_2_base p1 = cvit->source(), p2 = cvit->target();
		a += p1.x()*p2.y() - p2.x()*p1.y();
	}
	return 0.5*a;
}

CDT_2 mesh2Dpixel(unsigned nx, unsigned ny, double width, double height, const MesherData& meshParams)
{
	std::list<CDT_2::Point> seedPts; // Empty
	CDT_2 outMesh = mesh2D(getSegRep(nx, ny, width, height), seedPts, meshParams);
	for(CDT_2::All_faces_iterator it = outMesh.all_faces_begin(); it != outMesh.all_faces_end(); ++it)
		it->info().nesting_level = 0;
	return outMesh;
}

CDT_2 mesh2Dpixel(unsigned nx, unsigned ny, double width, double height, 
			const std::vector<GenericMaterial>& inMatVec, const MesherData& meshParams)
{
	std::list<CDT_2::Point> seedPts; // Empty
	CDT_2 outMesh = mesh2D(getSegRep(nx, ny, width, height), seedPts, meshParams);
	// Map material properties to triangles
	// First set all nesting_levels
	for(CDT_2::All_faces_iterator it = outMesh.all_faces_begin(); it != outMesh.all_faces_end(); ++it)
		it->info().nesting_level = 0;
	// Now set materials for finite faces
	double dx = width/(double)nx, dy = height/(double)ny;
	for(CDT_2::Finite_faces_iterator it = outMesh.finite_faces_begin(); it != outMesh.finite_faces_end(); ++it)
	{
		it->info().nesting_level = 1;
		// Get centroid of triangle and determine which pixel it falls in
		Point_2_base p1 = it->vertex(0)->point();
		Point_2_base p2 = it->vertex(1)->point();
		Point_2_base p3 = it->vertex(2)->point();
		Point_2_base pc = CGAL::centroid(p1, p2, p3);
		// TODO: For now this is assumed to start at (0,0)
		unsigned kx = (unsigned)(pc.x()/dx), ky = (unsigned)(pc.y()/dy);
		it->info().triMat = inMatVec[ky + kx*ny];
		it->info().optVal = 1.;
	}
	return outMesh;
}

CDT_2 mesh2Dpixel(unsigned nx, unsigned ny, double width, double height,
									const std::vector<double>& inOptValVec, const MesherData& meshParams)
{
	std::list<CDT_2::Point> seedPts; // Empty
	CDT_2 outMesh = mesh2D(getSegRep(nx, ny, width, height), seedPts, meshParams);
	// Map material properties to triangles
	// First set all nesting_levels
	for(CDT_2::All_faces_iterator it = outMesh.all_faces_begin(); it != outMesh.all_faces_end(); ++it)
		it->info().nesting_level = 0;
	// Now set materials for finite faces
	double dx = width/(double)nx, dy = height/(double)ny;
	for(CDT_2::Finite_faces_iterator it = outMesh.finite_faces_begin(); it != outMesh.finite_faces_end(); ++it)
	{
		it->info().nesting_level = 1;
		// Get centroid of triangle and determine which pixel it falls in
		Point_2_base p1 = it->vertex(0)->point();
		Point_2_base p2 = it->vertex(1)->point();
		Point_2_base p3 = it->vertex(2)->point();
		Point_2_base pc = CGAL::centroid(p1, p2, p3);
		// TODO: For now this is assumed to start at (0,0)
		unsigned kx = (unsigned)(pc.x()/dx), ky = (unsigned)(pc.y()/dy);
		it->info().optVal = inOptValVec[ky + kx*ny];
	}
	return outMesh;
}

CDT_2 mesh2D(const std::vector<Mesh_Segment_2>& inSegVec, const MesherData& meshParams)
{
	std::vector<GenericMaterial> emptyMatVec;
	return mesh2D(inSegVec, emptyMatVec, meshParams);
}

CDT_2 mesh2D(const std::vector<std::vector<Mesh_Segment_2> >& inSegVec, const MesherData& meshParams)
{
	std::vector<GenericMaterial> emptyMatVec;
	return mesh2D(inSegVec, emptyMatVec, meshParams);
}

CDT_2 mesh2D(const std::vector<Mesh_Segment_2>& inSegVec, const std::vector<GenericMaterial>& inMatVec,
              const MesherData& meshParams)
{
	std::list<CDT_2::Point> seedPts; // Empty
	CDT_2 outMesh = mesh2D(inSegVec, seedPts, meshParams);
	// Map material properties to triangles
	mark_domains(outMesh);
	if(inMatVec.size() == 1)
		setTriangleInfo(outMesh, inMatVec[0]);
	else if(inMatVec.size() > 1)
		setTriangleInfo(outMesh, inMatVec[0], inMatVec[1]);
	return outMesh;
}

CDT_2 mesh2D(const std::vector<std::vector<Mesh_Segment_2> >& inSegVec, const std::vector<GenericMaterial>& inMatVec,
              const MesherData& meshParams)
{
	std::list<CDT_2::Point> seedPts;
	computeSeedPoints(inSegVec, seedPts);
	std::vector<Mesh_Segment_2> flatSegVec;
	flattenMeshSegments(inSegVec, flatSegVec);
	CDT_2 outMesh = mesh2D(flatSegVec, seedPts, meshParams);
	// Map material properties to triangles
	mark_domains(outMesh);
	if(inMatVec.size() == 1)
		setTriangleInfo(outMesh, inMatVec[0]);
	else if(inMatVec.size() > 1)
		setTriangleInfo(outMesh, inMatVec[0], inMatVec[1]);
	return outMesh;
}

C3t3 mesh3D(const Mesh_polyhedron_3& polyToMesh, const MesherData& mp)
{
	C3t3 outMesh;
	// Set mesh criteria
	Mesh_criteria criteria(mp.tetMeshEdgeSize, mp.tetMeshFacetAngle, 
												 mp.tetMeshFacetSize, mp.tetMeshFacetDistance,
                         mp.tetMeshCellRadiusEdgeRatio, mp.tetMeshCellSize);
	// Create domain
	Mesh_domain domain(polyToMesh);
	domain.detect_features(10.); // Ensures that sharp edges are preserved
	// Generate mesh
	outMesh = CGAL::make_mesh_3<C3t3>(domain, criteria);
	return outMesh;
}

PolyhedronBuilderFromTOMesh::PolyhedronBuilderFromTOMesh(TOMesh3DSurface tom):
	upTOM(std::unique_ptr<TOMesh3DSurface>(new TOMesh3DSurface(std::move(tom))))
{
}

void PolyhedronBuilderFromTOMesh::operator()(HalfedgeDS& hds)
{
	CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(hds, true);
	B.begin_surface(upTOM->getNumNodes(), upTOM->getNumElements());
	// Copy points
	for(std::size_t k = 0; k < upTOM->getNumNodes(); ++k)
		B.add_vertex(upTOM->getNode3D(k));
	// Copy facets
	for(std::size_t kf = 0; kf < upTOM->getNumElements(); ++kf)
	{
		B.begin_facet();
		const std::vector<std::size_t>& curElem = upTOM->getElementConnectivity(kf);
		for(auto kit = curElem.begin(); kit != curElem.end(); ++kit)
			B.add_vertex_to_facet(*kit);
		B.end_facet();
	}
	if(B.error())
		B.rollback();
	B.end_surface();
}

void addPointsToSegments(const std::vector<Point_2_base>& ptVec, std::vector<Mesh_Segment_2>& segVec)
{
	std::vector<Mesh_Segment_2> tmpVec;
	generateSegment2FromPoints(ptVec, tmpVec);
	segVec.insert(segVec.end(), tmpVec.begin(), tmpVec.end());
}

void addPointsToSegments(const std::vector<Point_2_base>& ptVec, std::vector<std::vector<Mesh_Segment_2> >& segVV)
{
	std::vector<Mesh_Segment_2> segVec;
	generateSegment2FromPoints(ptVec, segVec);
	segVV.push_back(segVec);
}

bool checkValidityAsBoundaryMesh(const std::vector<Mesh_Segment_2>& segVec)
{
	CGAL::Polygon_2<Mesh_K> polygon;
	meshSegs2Polygon(segVec, polygon);
	return polygon.is_simple();
}

void generateSegment2FromPoints(const std::vector<Point_2_base>& ptVec, std::vector<Mesh_Segment_2>& segVec)
{
	segVec.resize(ptVec.size());
	for(unsigned k = 0; k < ptVec.size(); ++k)
	{
		unsigned kp1 = (k + 1) % ptVec.size();
		Point_2_base p1 = ptVec[k], p2 = ptVec[kp1];
		Mesh_Segment_2 tmpSeg(p1, p2);
		segVec[k] = tmpSeg;
	}
}

void generatePointsFromSegment2(const std::vector<Mesh_Segment_2>& segVec, std::vector<Point_2_base>& ptVec)
{
	ptVec.resize(segVec.size());
	for(unsigned k = 0; k < segVec.size(); ++k)
		ptVec[k] = segVec[k].source();
}

void flattenMeshSegments(const std::vector<std::vector<Mesh_Segment_2> >& orderedVec, std::vector<Mesh_Segment_2>& outSegVec)
{
	outSegVec.clear();
	for(unsigned kshp = 0; kshp < orderedVec.size(); ++kshp)
		for(unsigned kseg = 0; kseg < orderedVec[kshp].size(); ++kseg)
			outSegVec.push_back(orderedVec[kshp][kseg]);
}

void orderMeshSegments(const std::vector<Mesh_Segment_2>& segVec, std::vector<std::vector<Mesh_Segment_2> >& orderedVecVec)
{
	vector<Mesh_Segment_2> unorderedVec = segVec;
	orderedVecVec.clear();
	if(segVec.size() != 0)
	{
		vector<Mesh_Segment_2> curRowVec;
		orderedVecVec.push_back(curRowVec);
		Uint curRow = 0;
		orderedVecVec[curRow].push_back(unorderedVec.back());
		unorderedVec.pop_back();
		while(!unorderedVec.empty())
		{
			vector<Mesh_Segment_2>::iterator vit;
			bool found = false, flip = false;
			Point_2_base searchPt = orderedVecVec[curRow].back().target();
			for(vit = unorderedVec.begin(); vit != unorderedVec.end() && !found; ++vit)
			{
				if(TOL_EQ(searchPt.x(), vit->source().x(), 1e-8) && TOL_EQ(searchPt.y(), vit->source().y(), 1e-8))
					found = true;
				else if(TOL_EQ(searchPt.x(), vit->target().x(), 1e-8) && TOL_EQ(searchPt.y(), vit->target().y(), 1e-8))
				{
					found = true;
					flip = true;
				}
			}
			if(found)
			{
				vit--;
				Mesh_Segment_2 orderedSeg = *vit;
				if(flip)
					orderedSeg = orderedSeg.opposite();
				orderedVecVec[curRow].push_back(orderedSeg);
				unorderedVec.erase(vit);
			}
			else
			{
				vector<Mesh_Segment_2> newRow;
				newRow.push_back(unorderedVec.back());
				orderedVecVec.push_back(newRow);
				curRow++;
				unorderedVec.pop_back();
			}
		}
		for(Uint k = 0; k < orderedVecVec.size(); ++k)
		{
			Real a = computeSignedArea(orderedVecVec[k]);
			if(a < 0.)
			{
				vector<Mesh_Segment_2>::iterator vit;
				for(vit = orderedVecVec[k].begin(); vit != orderedVecVec[k].end(); ++vit)
					*vit = vit->opposite();
				std::reverse(orderedVecVec[k].begin(), orderedVecVec[k].end());
			}
		}
	}
}

TOMesh3DSurface linearExtrude(const CDT_2& inMesh, double length)
{
	// Copy vertices to TOMesh data format
	std::vector<Point_3_base> ptVec(2*inMesh.number_of_vertices());
	std::map<CDT_2::Vertex_handle, std::size_t> vitMap;
	std::size_t kp = 0;
	for(CDT_2::Finite_vertices_iterator vit = inMesh.finite_vertices_begin();
			vit != inMesh.finite_vertices_end(); ++vit, ++kp)
	{
		Point_2_base tmpp = vit->point();
		ptVec[kp] = Point_3_base(tmpp.x(), tmpp.y(), 0.);
		vitMap[vit] = kp;
	}
	for(CDT_2::Finite_vertices_iterator vit = inMesh.finite_vertices_begin();
			vit != inMesh.finite_vertices_end(); ++vit, ++kp)
	{
		Point_2_base tmpp = vit->point();
		ptVec[kp] = Point_3_base(tmpp.x(), tmpp.y(), length);
	}
	// Construct all facets
	std::vector<std::vector<std::size_t>> connectivityVec;
	connectivityVec.reserve(3*inMesh.number_of_faces());
	std::size_t offset = inMesh.number_of_vertices();
	for(CDT_2::Finite_faces_iterator fit=inMesh.finite_faces_begin();
			fit!=inMesh.finite_faces_end(); ++fit)
	{
		if(fit->info().in_domain())
		{
			connectivityVec.push_back({vitMap[fit->vertex(0)], vitMap[fit->vertex(1)], vitMap[fit->vertex(2)]});
			connectivityVec.push_back({vitMap[fit->vertex(0)]+offset, vitMap[fit->vertex(2)]+offset, vitMap[fit->vertex(1)]+offset});
		}
	}
	// Side facets, 2 triangles form the rectangles that connect the upper and lower planes
	for(CDT_2::Finite_edges_iterator eit = inMesh.finite_edges_begin();
			eit != inMesh.finite_edges_end(); ++eit)
	{
		if(inMesh.is_constrained(*eit))
		{
			CDT_2::Vertex_handle fVertex = eit->first->vertex(CDT_2::ccw(eit->second));
			CDT_2::Vertex_handle sVertex = eit->first->vertex(CDT_2::cw(eit->second));
			connectivityVec.push_back({vitMap[fVertex], vitMap[sVertex], vitMap[sVertex] + offset});
			connectivityVec.push_back({vitMap[sVertex] + offset, vitMap[fVertex] + offset, vitMap[fVertex]});
		}
	}
	return TOMesh3DSurface(ptVec, connectivityVec);
}

}// namespace GeometryTranslation
}// namespace Topologies
