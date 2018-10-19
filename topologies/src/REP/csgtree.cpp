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

#include "csgtree.h"
#include "geometrytranslation.h"
#include "helper.h"
#include "filter2d.h"
#include "tomesh.h"
#include "csgnode.h"
#include "csgfunctionnode.h"
#include "csgterminalnode.h"
#include <algorithm>
#include <functional>
#include <CGAL/Nef_2/Object_handle.h>

namespace Topologies{
const double CSGTreeRep::maxScale = 3., CSGTreeRep::maxMove = 0.5;

namespace
{
	// Nonmember helper funcs
	double nefPolyArea(Nef_polyhedron_2& inNef)
	{
		Explorer_2 E = inNef.explorer();
		Explorer_2::Face_const_iterator f = E.faces_begin();
		double nefa = 0.;
		for(++f; f != E.faces_end(); ++f)
		{
			Explorer_2::Halfedge_around_face_const_circulator he_circ = E.face_cycle(f);
			if(!E.is_standard(he_circ->vertex()))
				continue;
			bool includeFace = f->mark();
			double facea = 0.;
			do
			{
				Point_2 tmp = E.point(he_circ->vertex()), tmp2 = E.point((++he_circ)->vertex());
				he_circ--;
				facea += tmp[0].to_double()*tmp2[1].to_double() - tmp[1].to_double()*tmp2[0].to_double();
				++he_circ;
			} while(he_circ != E.face_cycle(f));
			if(includeFace)
				nefa += fabs(facea);
			else
				nefa -= fabs(facea);
		}
		return 0.5*nefa;
	}

	bool isPointInNefPoly(const Nef_polyhedron_2& nefPoly, const Point_2& p)
	{
		CGAL::Object_handle h = nefPoly.locate(p);
		return nefPoly.contains(h);
	}

	double getIntersectionPercent(const Nef_polyhedron_2& nefPoly, CDT_2::Finite_faces_iterator& fit)
	{
		// Construct new Nef_polyhedron_2 out of face iterator
		std::vector<Point_2> ptVec(3);
		for(unsigned k = 0; k < 3; ++k)
			ptVec[k] = Point_2(fit->vertex(k)->point().x(), fit->vertex(k)->point().y());
		Nef_polyhedron_2 Ntri(ptVec.begin(), ptVec.end());
		Nef_polyhedron_2 Nint = nefPoly.intersection(Ntri);
		double a1 = nefPolyArea(Nint);
		FT a2 = CGAL::area(ptVec[0], ptVec[1], ptVec[2]);
		return a1/a2.to_double();
	}
	
	std::list<W_Point_2> vec2list(const std::vector<Mesh_Segment_2>& inSegVec)
	{
		std::vector<Point_2_base> ptVec;
		GeometryTranslation::generatePointsFromSegment2(inSegVec, ptVec);
		std::list<W_Point_2> outList;
		for(auto it = ptVec.begin(); it != ptVec.end(); ++it)
			outList.push_back(W_Point_2(*it, 0.));
		return outList;
	}
}

CSGTreeRep::CSGTreeRep(const InputLoader::TORCSGTree& inputParams) :
	TopOptRep(tortCSG2D),
	minDensity(inputParams.getMinDensity()),
	width(inputParams.getRegionDimensions(0)),
	height(inputParams.getRegionDimensions(1)),
	nx(inputParams.getMeshSizes(0)),
	ny(inputParams.getMeshSizes(1)),
	nShapesX(inputParams.getShapeNums(0)),
	nShapesY(inputParams.getShapeNums(1)),
	nPointsPerShape(inputParams.getNumPointsPerShape()),
	useAffine(inputParams.getUseAffine()),
	shapesAreHoles(inputParams.getShapesAreHoles())
{
	finishSetup();
	
}

CSGTreeRep::CSGTreeRep(const std::vector<std::vector<int> >& discreteParams, const std::vector<std::vector<double> >& realParams) :
	TopOptRep(tortCSG2D)
{
	bool error = discreteParams.empty();
	if(!error)
	{
		error = discreteParams[0].size() != 7;
		if(!error)
		{
			nx = discreteParams[0][0];
			ny = discreteParams[0][1];
			nShapesX = discreteParams[0][2];
			nShapesY = discreteParams[0][3];
			nPointsPerShape = discreteParams[0][4];
			useAffine = (bool)discreteParams[0][5];
			shapesAreHoles = (bool)discreteParams[0][6];
		}
	}
	error |= realParams.empty();
	if(!error)
	{
		error = realParams[0].size() != 3;
		if(!error)
		{
			minDensity = realParams[0][0];
			width = realParams[0][1];
			height = realParams[0][2];
		}
	}
	if(error)
	{
		std::cout << "Error in constructor CSGTreeRep, constructing from vectors of defining parameters." << std::endl;
		std::cout << "Constructor used for MPI code, check there for error!" << std::endl;
		abort();
	}
	finishSetup();
}

CSGTreeRep::CSGTreeRep(const InputLoader::TORCSGTree& inputParams, const std::vector<Mesh_Segment_2>& segVec) :
	TopOptRep(tortCSG2D),
	minDensity(inputParams.getMinDensity()),
	width(inputParams.getRegionDimensions(0)),
	height(inputParams.getRegionDimensions(1)),
	nx(inputParams.getMeshSizes(0)),
	ny(inputParams.getMeshSizes(1)),
	nShapesX(inputParams.getShapeNums(0)),
	nShapesY(inputParams.getShapeNums(1)),
	nPointsPerShape(inputParams.getNumPointsPerShape()),
	useAffine(inputParams.getUseAffine()),
	shapesAreHoles(inputParams.getShapesAreHoles())
{
	finishSetup(segVec);
}

void CSGTreeRep::finishSetup()
{
	finishSetup({});
}

void CSGTreeRep::finishSetup(const std::vector<Mesh_Segment_2>& segVec)
{
	std::vector<Point_2> boundPoly(4);
	boundPoly[0] = Point_2(0., 0.);
	boundPoly[1] = Point_2(width, 0.);
	boundPoly[2] = Point_2(width, height);
	boundPoly[3] = Point_2(0., height);
	CSGNode::setDecodeType(dtBoolean);
	CSGNode::setDimensionType(dt2d);
	CSGTerminalNode::setBoundPoly(boundPoly);
	CSGTerminalNode::setUseAlpha(false);
	if(segVec.empty())
	{
		if(nShapesY*nShapesX == 1)
			rootNode = std::unique_ptr<CSGNode>(new CSGTerminalNode(0, nPointsPerShape, Point_2_base(width*0.5, height*0.5),
																				0.25*width, 0.25*height));
		else
			rootNode = std::unique_ptr<CSGNode>(new CSGFunctionNode(0, nShapesX, nShapesY, nPointsPerShape, width, height));
	}
	else
	{
		std::vector<std::vector<Mesh_Segment_2>> segVV;
		GeometryTranslation::partitionMeshSegPolys(segVec, segVV);
		std::vector<std::list<W_Point_2>> ptListVec(segVV.size());
		for(std::size_t k = 0; k < segVV.size(); ++k)
			ptListVec[k] = vec2list(segVV[k]);
		rootNode = std::unique_ptr<CSGNode>(new CSGFunctionNode(0, ptListVec));
	}
	mesh = GeometryTranslation::mesh2Dpixel(nx, ny, width, height, getDefaultMeshParams());
	if(useAffine)
	{
		std::vector<CSGNode*> nodeVec;
		rootNode->getTerminalNodeList(nodeVec);
		affineVec.resize(3*nodeVec.size());
		for(std::size_t k = 0; k < nodeVec.size(); ++k)
		{
			affineVec[3*k    ] = 1.; // Scale
			affineVec[3*k + 1] = 0.; // X translate
			affineVec[3*k + 2] = 0.; // Y translate
		}
	}
}

CSGTreeRep::CSGTreeRep(const CSGTreeRep& copy) :
	TopOptRep(copy),
	minDensity(copy.minDensity),
	width(copy.width),
	height(copy.height),
	nx(copy.nx),
	ny(copy.ny),
	nShapesX(copy.nShapesX),
	nShapesY(copy.nShapesY),
	nPointsPerShape(copy.nPointsPerShape),
	mesh(copy.mesh),
	useAffine(copy.useAffine),
	affineVec(copy.affineVec),
	shapesAreHoles(copy.shapesAreHoles)
{
	if(copy.rootNode->getNumChildren() == 0)
		rootNode = std::unique_ptr<CSGNode>(new CSGTerminalNode(*dynamic_cast<CSGTerminalNode*>(copy.rootNode.get())));
	else
		rootNode = std::unique_ptr<CSGNode>(new CSGFunctionNode(*dynamic_cast<CSGFunctionNode*>(copy.rootNode.get())));
}

CSGTreeRep::CSGTreeRep(CSGTreeRep && copy) : 
	TopOptRep(copy) 
{
	swap(copy);
}

CSGTreeRep& CSGTreeRep::operator=(CSGTreeRep copy)
{
	swap(copy); 
	return *this;
}

CSGTreeRep::~CSGTreeRep()
{
}

void CSGTreeRep::swap(CSGTreeRep& arg2)
{
	TopOptRep::swap(arg2);
	std::swap(minDensity, arg2.minDensity);
	std::swap(width, arg2.width);
	std::swap(height, arg2.height);
	std::swap(nx, arg2.nx);
	std::swap(ny, arg2.ny);
	std::swap(nPointsPerShape, arg2.nPointsPerShape);
	std::swap(nShapesX, arg2.nShapesX);
	std::swap(nShapesY, arg2.nShapesY);
	mesh.swap(arg2.mesh);
	std::swap(rootNode, arg2.rootNode);
	std::swap(useAffine, arg2.useAffine);
	affineVec.swap(arg2.affineVec);
	std::swap(shapesAreHoles, arg2.shapesAreHoles);
}

std::unique_ptr<TopOptRep> CSGTreeRep::clone() const
{
	std::unique_ptr<TopOptRep> outUP(new CSGTreeRep(*this));
	return outUP;
}

void CSGTreeRep::initialize()
{
}	

void CSGTreeRep::initialize(double val, std::pair<double, double> randRange)
{
	std::vector<CSGNode*> nodeVec;
	rootNode->getTerminalNodeList(nodeVec);
	// Apply random translations and scaling to points:
	for(std::size_t k = 0, kp = 0; k < nodeVec.size(); ++k)
	{
		CSGTerminalNode* curNode = dynamic_cast<CSGTerminalNode*>(nodeVec[k]);
		double scaleFactor = HelperNS::RandomGen::instance().randRealInRange<double>(0.1, 0.5);
		curNode->scalePoints(scaleFactor*maxScale);
		double xmv = HelperNS::RandomGen::instance().randRealInRange(randRange.first, randRange.second),
					ymv = HelperNS::RandomGen::instance().randRealInRange(randRange.first, randRange.second);
		Vector_2 mvv2(xmv, ymv);
		curNode->translatePoints(mvv2);
	}
}

void CSGTreeRep::initialize(double val)
{
}

void CSGTreeRep::randomize()
{
	if(!useAffine)
	{
		std::vector<CSGNode*> nodeVec;
		rootNode->getTerminalNodeList(nodeVec);
		// Apply random translations and scaling to points:
		for(std::size_t k = 0, kp = 0; k < nodeVec.size(); ++k)
		{
			CSGTerminalNode* curNode = dynamic_cast<CSGTerminalNode*>(nodeVec[k]);
			double scaleFactor = HelperNS::RandomGen::instance().randRealInRange(1./maxScale, maxScale);
			curNode->scalePoints(scaleFactor);
			double xmv = HelperNS::RandomGen::instance().randRealInRange(-maxMove*width, maxMove*width),
						 ymv = HelperNS::RandomGen::instance().randRealInRange(-maxMove*height, maxMove*height);
			Vector_2 mvv2(xmv, ymv);
			curNode->translatePoints(mvv2);
		}
	}
	else
	{
		std::vector<CSGNode*> nodeVec;
		rootNode->getTerminalNodeList(nodeVec);
		for(std::size_t k = 0, kp = 0; k < nodeVec.size(); ++k)
		{
			double &scale = affineVec[3*k], &xt = affineVec[3*k + 1], &yt = affineVec[3*k + 2];
			scale = HelperNS::RandomGen::instance().randRealInRange(1./maxScale, maxScale);
			xt = HelperNS::RandomGen::instance().randRealInRange(-maxMove*width, maxMove*width);
			yt = HelperNS::RandomGen::instance().randRealInRange(-maxMove*width, maxMove*width);
		}
	}
}

// Decode
std::unique_ptr<TOMesh> CSGTreeRep::get3DSurfaceMesh() const {return nullptr;}
std::unique_ptr<TOMesh> CSGTreeRep::get3DVolumeMesh() const {return nullptr;}

void CSGTreeRep::get2DSegments(std::vector<Mesh_Segment_2>& segVec) const
{
	Nef_polyhedron_2 nefp = getNefPolyWithBoundaries();
	Explorer_2 E = nefp.explorer();
	Explorer_2::Face_const_iterator f = E.faces_begin();
	for(++f; f != E.faces_end(); ++f)
	{
		Explorer_2::Halfedge_around_face_const_circulator he_circ = E.face_cycle(f);
		if(!E.is_standard(he_circ->vertex()))
			continue;
		for(bool first = true; he_circ != E.face_cycle(f) || first; ++he_circ, first = false)
		{
			Point_2 p1 = E.point(he_circ->vertex());
			--he_circ;
			Point_2 p2 = E.point(he_circ->vertex());
			++he_circ;
			Point_2_base pb1(p1.x().to_double(), p1.y().to_double());
			Point_2_base pb2(p2.x().to_double(), p2.y().to_double());
			Mesh_Segment_2 seg(pb1, pb2);
			segVec.push_back(seg);
		}
	}
//	std::vector<Mesh_Segment_2> boundaryVec;
//	getBoundary(boundaryVec);
//	segVec.insert(segVec.end(), boundaryVec.begin(), boundaryVec.end());
}

GeometryTranslation::MesherData CSGTreeRep::getDefaultMeshParams() const
{
	GeometryTranslation::MesherData meshParams;
	meshParams.triMeshEdgeAngle = 10.;
	meshParams.triMeshEdgeSize = width/(double)nx;
	return meshParams;
}

std::unique_ptr<TOMesh> CSGTreeRep::get2DMesh() const
{
	return get2DMesh(getDefaultMeshParams());
}

std::unique_ptr<TOMesh> CSGTreeRep::get2DMesh(const GeometryTranslation::MesherData& meshParams) const
{
	setMeshOptVals();
	return std::unique_ptr<TOMesh>(new TOMesh2D(mesh));
}

std::unique_ptr<TOMesh> CSGTreeRep::getOutputMesh() const
{
	setMeshOptVals();
	return std::unique_ptr<TOMesh>(new TOMesh2D(mesh));
}

void CSGTreeRep::getBoundary(std::vector<Mesh_Segment_2>& boundaryVec) const
{
	boundaryVec.resize(4);
	boundaryVec[0] = Mesh_Segment_2(Point_2_base(0., 0.), Point_2_base(width, 0.));
	boundaryVec[1] = Mesh_Segment_2(Point_2_base(width, 0.), Point_2_base(width, height));
	boundaryVec[2] = Mesh_Segment_2(Point_2_base(width, height), Point_2_base(0., height));
	boundaryVec[3] = Mesh_Segment_2(Point_2_base(0., height), Point_2_base(0., 0.));
}

// Modify
void CSGTreeRep::refine()
{
	// TODO: Add more children to the tree?
}

void CSGTreeRep::prune()
{
	// TODO: Not implemented yet
}

void CSGTreeRep::scaleParameters(std::vector<double>& newvals) const
{
	if(!useAffine)
	{
		// Optimization parameters come in as between 0-1, scale them appropriately
		for(std::size_t k = 0; k < newvals.size(); k+=2)
		{
			newvals[k] = newvals[k]*width;
			newvals[k+1] = newvals[k+1]*height;
		}
	}
	else
	{
		for(std::size_t k = 0; k < newvals.size(); k+=3)
		{
			newvals[k] = newvals[k]*maxScale;
			newvals[k+1] = maxMove*width*(2.*newvals[k+1] - 1.);
			newvals[k+2] = maxMove*height*(2.*newvals[k+2] - 1.);
		}
	}
}

void CSGTreeRep::unscaleParameters(std::vector<double>& newvals) const
{
	if(!useAffine)
	{
		// Optimization parameters come in as between 0-1, scale them appropriately
		for(std::size_t k = 0; k < newvals.size(); k+=2)
		{
			newvals[k] = newvals[k]/width;
			newvals[k+1] = newvals[k+1]/height;
		}
	}
	else
	{
		for(std::size_t k = 0; k < newvals.size(); k+=3)
		{
			newvals[k] = newvals[k]/maxScale;
			newvals[k+1] = 0.5*(newvals[k+1]/width/maxMove + 1.);
			newvals[k+2] = 0.5*(newvals[k+2]/height/maxMove + 1.);
		}
	}
}

// Data access
void CSGTreeRep::updateRealRep()
{
	affineVec = realOptVals;
	scaleParameters(affineVec);
	if(!useAffine)
	{
		std::vector<double>::const_iterator vit = affineVec.begin();
		rootNode->setVector(vit, affineVec.end());
	}
}

void CSGTreeRep::setDiscreteRep(const std::vector<int>& newvals)
{
}

void CSGTreeRep::setMPIRep(const std::vector<std::vector<int>>& discreteVars, const std::vector<std::vector<double>>& realVars)
{
	// First check realVars size.  
	// The MPI code may send only the realRep rather than the full MPIRep
	// This is a little clunky now, and should be fixed
	if(realVars.size() == 1)
		setRealRep(realVars[0].begin(), realVars[0].end());
	else
	{
		assert(discreteVars.size() == 2);
		assert(realVars.size() == 4);
		std::vector<double>::const_iterator ptIt = realVars[0].begin(), matIt = realVars[1].begin(), prIt = realVars[2].begin();
		std::vector<int>::const_iterator treeIt = discreteVars[0].begin(), nmIt = discreteVars[1].begin();
		int nodeType = *treeIt;
		++treeIt;
		if(nodeType >= 0)
			rootNode = std::unique_ptr<CSGNode>(new CSGFunctionNode(0, nodeType, nmIt, ptIt, matIt, prIt, treeIt));
		else
			rootNode = std::unique_ptr<CSGNode>(new CSGTerminalNode(0, *nmIt, ptIt, matIt, prIt, treeIt));
		affineVec = realVars[3];
	}
}

void CSGTreeRep::getRealRep(std::vector<double>& realVec) const
{
	if(!useAffine)
	{
		std::list<W_Point_2> ptList;
		rootNode->getLocalOptFormat(ptList);
		realVec.resize(2*ptList.size());
		std::size_t k = 0;
		for(std::list<W_Point_2>::const_iterator lit = ptList.begin(); lit != ptList.end(); ++lit)
		{
			realVec[k++] = lit->x();
			realVec[k++] = lit->y();
		}
	}
	else
		realVec = affineVec;
//	std::cout << "unscaling-----------------------------" << std::endl;
//	std::cout << "vals_beforeUnscale = [";
//  for(auto it = realVec.begin(); it != realVec.end(); ++it)
//    std::cout << *it << " ";
//  std::cout << "];" << std::endl;
	unscaleParameters(realVec); // Set all parameters to 0-1 range
//	std::cout << "vals_after = [";
//  for(auto it = realVec.begin(); it != realVec.end(); ++it)
//    std::cout << *it << " ";
//  std::cout << "];" << std::endl;
}

void CSGTreeRep::getDiscreteRep(std::vector<int>& discVec) const
{
}

void CSGTreeRep::getMPIRep(std::vector<std::vector<int> >& discreteVars, std::vector<std::vector<double> >& realVars) const
{
	std::vector<int> treeVec;
	std::vector<W_Point_2> ptVec;
	std::vector<GenericMaterial> matVec;
	std::vector<Real> priorityVec;
	rootNode->getMPIDataFormat(treeVec, ptVec, matVec, priorityVec);
	std::vector<double> flatPtVec(ptVec.size()*3);
	for(std::size_t k = 0; k < ptVec.size(); ++k)
	{
		flatPtVec[3*k] = ptVec[k].x();
		flatPtVec[3*k + 1] = ptVec[k].y();
		flatPtVec[3*k + 2] = ptVec[k].weight();
	}
	std::vector<double> flatMatVec;
	std::vector<int> numMatParamVec(matVec.size());
	for(std::size_t k = 0; k < matVec.size(); ++k)
	{
		unsigned nm = matVec[k].getNumParameters();
		numMatParamVec[k] = nm;
		for(unsigned km = 0; km < nm; ++km)
			flatMatVec.push_back(matVec[k].getParameter(km));
	}
	realVars.resize(4);
	realVars[0] = flatPtVec;
	realVars[1] = flatMatVec;
	realVars[2] = priorityVec;
	realVars[3] = affineVec;
	discreteVars.resize(2);
	discreteVars[0] = treeVec;
	discreteVars[1] = numMatParamVec;
}

std::size_t CSGTreeRep::getDataSize() const
{
	if(!useAffine)
		return rootNode->countPoints()*2;
	else
		return affineVec.size();
}

void CSGTreeRep::getDataSize(std::vector<std::size_t>& sizes) const
{
	sizes.resize(2);
	sizes[0] = getDataSize();
	sizes[1] = 1;
}

double CSGTreeRep::computeVolumeFraction() const
{
	Nef_polyhedron_2 Nint = getNefPolyWithBoundaries();
	double a1 = nefPolyArea(Nint);
	if(shapesAreHoles)
		a1 = width*height - a1;
	double a2 = width*height;
	return a1/a2;
}

void CSGTreeRep::getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
                                     std::vector<std::vector<double> >& realParams) const
{
	discreteParams.resize(1);
	discreteParams[0] = {(int)nx, (int)ny, (int)nShapesX, (int)nShapesY, (int)nPointsPerShape, (int)useAffine, (int)shapesAreHoles};
	realParams.resize(1);
	realParams[0] = {minDensity, width, height};
}

void CSGTreeRep::filterData(std::vector<double>& valVec, double radius) const
{
}

void CSGTreeRep::filterData(double radius)
{
}

Nef_polyhedron_2 CSGTreeRep::getNefPoly() const
{
	if(!useAffine)
	{
		Nef_polyhedron_2 outNP;
		rootNode->decode(outNP);
		return outNP;
	}
	// use affine transformations
	std::unique_ptr<CSGNode> rootCopy;
	if(rootNode->getNumChildren() == 0)
		rootCopy = std::unique_ptr<CSGNode>(new CSGTerminalNode(*dynamic_cast<CSGTerminalNode*>(rootNode.get())));
	else
		rootCopy = std::unique_ptr<CSGNode>(new CSGFunctionNode(*dynamic_cast<CSGFunctionNode*>(rootNode.get())));
	std::vector<CSGNode*> nodeVec;
	rootCopy->getTerminalNodeList(nodeVec);
	for(std::size_t k = 0; k < nodeVec.size(); ++k)
	{
		double scale = affineVec[3*k], xt = affineVec[3*k + 1], yt = affineVec[3*k + 2];
		dynamic_cast<CSGTerminalNode*>(nodeVec[k])->scalePoints(scale);
		dynamic_cast<CSGTerminalNode*>(nodeVec[k])->translatePoints(Vector_2(xt, yt));
	}
	Nef_polyhedron_2 outNP;
	rootCopy->decode(outNP);
	return outNP;
}

Nef_polyhedron_2 CSGTreeRep::getNefPolyWithBoundaries() const
{
	Nef_polyhedron_2 nefp = getNefPoly();
	std::vector<Point_2> ptVec = {Point_2(0., 0.), Point_2(width, 0.), Point_2(width, height), Point_2(0., height)};
	Nef_polyhedron_2 Nbnds(ptVec.begin(), ptVec.end());
	return Nbnds.intersection(nefp);
}

void CSGTreeRep::setMeshOptVals() const
{
	Nef_polyhedron_2 nefPoly = getNefPoly();
	for(CDT_2::Finite_faces_iterator it = mesh.finite_faces_begin(); it != mesh.finite_faces_end(); ++it)
	{
		if(it->is_in_domain())
		{
			// First check if all vertices are in the CSG shape
			bool allin = true, nonein = true;
			for(unsigned k = 0; k < 3; ++k)
			{
				Point_2 p(it->vertex(k)->point().x(), it->vertex(k)->point().y());
				bool curPtin = isPointInNefPoly(nefPoly, p);
				allin &= curPtin;
				nonein &= !curPtin;
			}
			it->info().optVal = 1.;
			if(nonein)
				it->info().optVal = minDensity;
			if(!allin)
				it->info().optVal = getIntersectionPercent(nefPoly, it);
			if(it->info().optVal < minDensity)
				it->info().optVal = minDensity;
			if(shapesAreHoles)
				it->info().optVal = 1. - it->info().optVal;
		}
	}
}
}

