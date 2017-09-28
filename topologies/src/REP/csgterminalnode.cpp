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

#include "csgterminalnode.h"
#include "geometrytranslation.h"

namespace Topologies{
std::vector<Point_2> CSGTerminalNode::boundingPoly;
Real CSGTerminalNode::alphamax = 1., CSGTerminalNode::weightmax = 0.;
bool CSGTerminalNode::useAlpha = false;

using namespace std;

void CSGTerminalNode::initMaterials()
{
	GenericMaterial backgroundMat, firstMat;
	if(materialList.size() > 0)
	{
        backgroundMat = materialList[0];
        firstMat = materialList[0];
    }
	if(materialList.size() > 1)
		firstMat = materialList[1];

	if(theDT == dtBoolean)
		itsGenericMaterial = firstMat;
	else if(theDT == dtRandMat)
	{
		itsGenericMaterial = backgroundMat;
		itsGenericMaterial.genRandomMat();
	}	
	else
	{
		std::size_t rid = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, materialList.size() - 1);
		itsGenericMaterial = materialList[rid];
	}
}

CSGTerminalNode::CSGTerminalNode(int curTreeDepth):
  CSGNode(curTreeDepth)
{
	initMaterials();
	priority = HelperNS::RandomGen::instance().randRealInRange(2.0, prioritymax);
	alpha = HelperNS::RandomGen::instance().randRealInRange(0., alphamax);
	if(!useAlpha)
		setAlphaAndWeightsConvex();
}

CSGTerminalNode::CSGTerminalNode(int curTreeDepth, unsigned numPts, Point_2_base center, double rx, double ry) :
	CSGNode(curTreeDepth)
{
	initMaterials();
	// Set up points
	for(unsigned k = 0; k < numPts; ++k)
	{
		double th = (double)k/(double)numPts * 2.*PI;
		double x = rx*cos(th) + center.x(), y = ry*sin(th) + center.y();
		geneList.push_back(W_Point_2(Point_2_base(x, y), 0.));
	}
	priority = 1.;
	alpha = alphamax;
	if(!useAlpha)
		setAlphaAndWeightsConvex();
}

CSGTerminalNode::CSGTerminalNode(const CSGTerminalNode& copy):
	CSGNode(copy),
	itsGenericMaterial(copy.itsGenericMaterial),
	priority(copy.priority),
	alpha(copy.alpha)
{
	for(std::list<W_Point_2>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
		geneList.push_back(*plist);
}

CSGTerminalNode::CSGTerminalNode(CSGTerminalNode && copy) :
	CSGNode(copy)
{
	swap(copy);
}

CSGTerminalNode& CSGTerminalNode::operator=(CSGTerminalNode copy)
{
	swap(copy);
	return *this;
}

void CSGTerminalNode::swap(CSGTerminalNode& arg2)
{
	geneList.swap(arg2.geneList);
	std::swap(alpha, arg2.alpha);
	std::swap(priority, arg2.priority);
	itsGenericMaterial.swap(arg2.itsGenericMaterial);
}

CSGTerminalNode::CSGTerminalNode(const CSGTerminalNode& copy, MirrorType inMT):
  CSGNode(copy),
  itsGenericMaterial(copy.itsGenericMaterial),
  priority(copy.priority),
	alpha(copy.alpha)
{
	if(inMT == mtNone || inMT == mtZ)
	{
  	for(std::list<W_Point_2>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
    	geneList.push_back(*plist);
	}
	else if(inMT == mtX)
	{
		for(std::list<W_Point_2>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
			geneList.push_back(W_Point_2(Point_2_base(-plist->x(), plist->y()), plist->weight()));
	}
	else if(inMT == mtY)
	{
		for(std::list<W_Point_2>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
			geneList.push_back(W_Point_2(Point_2_base(plist->x(), -plist->y()),plist->weight()));
	}
}

CSGTerminalNode::CSGTerminalNode(int curTreeDepth, std::ifstream& saveFile):
	CSGNode(curTreeDepth)
{
	int numPoints, numMatParams;
	saveFile >> numMatParams;
	std::vector<Real> constParams;
	for(int k = 0; k < numMatParams; k++)
	{
		Real tmp;
		saveFile >> tmp;
		constParams.push_back(tmp);
	}
	std::vector<Real> rgb(3);
	saveFile >> rgb[0] >> rgb[1] >> rgb[2];
	GenericMaterial backgroundMat;
	if(materialList.size() > 0)
		backgroundMat = materialList[0];
	itsGenericMaterial = GenericMaterial(constParams, backgroundMat, rgb);
	saveFile >> priority;
	saveFile >> alpha;
	saveFile >> numPoints;
	for(int k = 0; k < numPoints; k++)
	{
		Real x, y, w;
		saveFile >> x >> y >> w;
		W_Point_2 newPoint(Point_2_base(x, y), w);
		boundsCheck(newPoint);
		geneList.push_back(newPoint);
	}
}

CSGTerminalNode::CSGTerminalNode(int curTreeDepth, int numMatParams, int& ptVecPos, int& matPos, int& prPos, int& treeVecPos, 
				Real* ptArray, Real* matArray, Real* prArray, int* treeArray):
	CSGNode(curTreeDepth)
{
	int numPoints = treeArray[treeVecPos++];
	std::vector<Real> constParams;
	for(int k = 0; k < numMatParams; k++)
		constParams.push_back(matArray[numMatParams*matPos + k]);
	matPos++;
	itsGenericMaterial = GenericMaterial(constParams);
	priority = prArray[prPos++];
	alpha = prArray[prPos++];
	for(int k = 0; k < numPoints; k++)
	{
		Real x, y, w;
		x = ptArray[ptVecPos++];
		y = ptArray[ptVecPos++];
		w = ptArray[ptVecPos++];
		W_Point_2 newPoint(Point_2_base(x, y), w);
		geneList.push_back(newPoint);
	}
}

CSGTerminalNode::CSGTerminalNode(int curTreeDepth, int numMatParams, std::vector<double>::const_iterator& ptIt, 
																std::vector<double>::const_iterator& matIt, std::vector<double>::const_iterator& prIt,
																std::vector<int>::const_iterator& treeIt) :
  CSGNode(curTreeDepth)
{
	std::vector<Real> constParams;
	for(int k = 0; k < numMatParams; k++)
	{
		constParams.push_back(*matIt);
		++matIt;
	}
	itsGenericMaterial = GenericMaterial(constParams);
	priority = *prIt;
	++prIt;
	alpha = *prIt;
	++prIt;
	unsigned numPoints = *treeIt;
	treeIt++;
	for(int k = 0; k < numPoints; k++)
	{
		Real x, y, w;
		x = *ptIt;
		++ptIt;
		y = *ptIt;
		++ptIt;
		w = *ptIt;
		++ptIt;
		W_Point_2 newPoint(Point_2_base(x, y), w);
		geneList.push_back(newPoint);
	}
}

CSGTerminalNode::CSGTerminalNode(int curTreeDepth, const std::list<W_Point_2>& ptList):
	CSGNode(curTreeDepth)
{
	geneList = ptList;
	GenericMaterial backgroundMat;
	if(materialList.size() > 0)
		backgroundMat = materialList[0];

	if(theDT == dtBoolean)
		itsGenericMaterial = backgroundMat;
	else if(theDT == dtRandMat)
	{
		itsGenericMaterial = backgroundMat;
		itsGenericMaterial.genRandomMat();
	}
	else
	{
		std::size_t rid = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, materialList.size()-1);
		itsGenericMaterial = materialList[rid];
	}
	priority = HelperNS::RandomGen::instance().randRealInRange(2.0, prioritymax);
	alpha = HelperNS::RandomGen::instance().randRealInRange(0., alphamax);
	if(!useAlpha)
		setAlphaAndWeightsConvex();
}

CSGTerminalNode::~CSGTerminalNode()
{
}

bool CSGTerminalNode::decodeMaterial(Point_2& testPt, GenericMaterial& theGM) const
{
	std::vector<Mesh_Segment_2> segVec;
	decode(segVec);
  return decodeMaterialFromSegs(testPt, segVec, theGM);
}

void CSGTerminalNode::decode(Nef_polyhedron_2& inNP) const
{
	std::vector<Mesh_Segment_2> segVec;
	decode(segVec);
	inNP = nefPolyFromSegs(segVec, Nef_polyhedron_2::EXCLUDED);
}

void CSGTerminalNode::decodeMask(Nef_polyhedron_2& inNP) const
{
	std::vector<Mesh_Segment_2> segVec;
	decode(segVec);
	inNP = nefPolyFromSegs(segVec, Nef_polyhedron_2::INCLUDED);
}

bool CSGTerminalNode::decodeMaterial(Point_3& testPt, GenericMaterial& theGM) const
{
	return false;
}

bool CSGTerminalNode::decodeMaterial(Point_3& testPt, int& theGM) const
{
	return false;
}

void CSGTerminalNode::setMaterialList(std::vector<GenericMaterial>& theGMVec) const
{
	theGMVec.push_back(itsGenericMaterial);
}

void CSGTerminalNode::decode(Nef_polyhedron_3& inNP) const
{
}

void CSGTerminalNode::decodeMask(Nef_polyhedron_3& inNP) const
{
}

void CSGTerminalNode::decode(std::vector<Mesh_Segment_2>& inSegVec) const
{
  std::list<W_Point_2> tmpGL;
	tmpGL = geneList;
	Alpha_shape_2 tmpA(tmpGL.begin(), tmpGL.end());
	tmpA.set_mode(Alpha_shape_2::REGULARIZED);
	tmpA.set_alpha(alpha);
	Alpha_shape_edges_iterator it;
	inSegVec.clear();
	for(it =  tmpA.alpha_shape_edges_begin(); it != tmpA.alpha_shape_edges_end(); ++it)
		inSegVec.push_back(tmpA.segment(*it));
}

void CSGTerminalNode::decode(std::vector< std::vector<Mesh_Segment_2> >& inSegVec) const
{
	std::list<W_Point_2> tmpGL;
	tmpGL = geneList;
	Alpha_shape_2 tmpA(tmpGL.begin(), tmpGL.end());
	tmpA.set_mode(Alpha_shape_2::REGULARIZED);
	inSegVec.clear();
	Alpha_shape_2::Alpha_iterator ait;
	for(ait = tmpA.alpha_begin(); ait != tmpA.alpha_end(); ++ait)
	{
		tmpA.set_alpha(*ait);
		Alpha_shape_edges_iterator it;
		vector<Mesh_Segment_2> tmpSegVec;
		for(it =  tmpA.alpha_shape_edges_begin(); it != tmpA.alpha_shape_edges_end(); ++it)
			tmpSegVec.push_back(tmpA.segment(*it));
		inSegVec.push_back(tmpSegVec);
	}
}

bool CSGTerminalNode::decodeMaterialFromSegs(Point_2& testPt, const std::vector<Mesh_Segment_2>& segVec,
                                                GenericMaterial& theGM) const
{
	std::vector<std::vector<Mesh_Segment_2> > orderedArray;
	GeometryTranslation::orderMeshSegments(segVec, orderedArray);
	for(Uint k = 0; k < orderedArray.size(); ++k)
	{
		vector<Point_2> ptPoly = getPointVecFromPoly(orderedArray[k]);
		if(CGAL::bounded_side_2(ptPoly.begin(), ptPoly.end(), testPt, K()) == CGAL::ON_BOUNDED_SIDE)
		{
			theGM = itsGenericMaterial;
			return true;
		}
	}
	return false;
}

Nef_polyhedron_2 CSGTerminalNode::nefPolyFromSegs(const std::vector<Mesh_Segment_2>& segVec,
                                                Nef_polyhedron_2::Boundary boundTreatment) const
{
	std::vector<std::vector<Mesh_Segment_2> > orderedArray;
	GeometryTranslation::orderMeshSegments(segVec, orderedArray);
	Nef_polyhedron_2 outPoly;
	for(Uint k = 0; k < orderedArray.size(); ++k)
	{
		vector<Point_2> ptPoly = getPointVecFromPoly(orderedArray[k]);
		Nef_polyhedron_2 curNP2(ptPoly.begin(), ptPoly.end(), boundTreatment);
		outPoly += curNP2;
	}
	return outPoly;
}

std::vector<Point_2> CSGTerminalNode::getPointVecFromPoly(const std::vector<Mesh_Segment_2>& orderedPoly) const
{
	// Copy all source points from Mesh_Segments
	vector<Point_2> outVec;
  for(Uint k = 0; k < orderedPoly.size(); ++k)
  {
    Mesh_K::Point_2 tmpMKP2 = orderedPoly[k].source();
    outVec.push_back(Point_2(tmpMKP2.x(), tmpMKP2.y()));
  }
	return outVec;
}

void CSGTerminalNode::movePointsOB()
{
	Real maxx = boundingPoly[0].x().to_double(), maxy = boundingPoly[0].y().to_double();
	for(unsigned int k1 = 0; k1 < boundingPoly.size(); k1++)
	{
		maxx = MAX(maxx, boundingPoly[k1].x().to_double());
		maxy = MAX(maxy, boundingPoly[k1].y().to_double());
	}
	Vector_2 tv(10*maxx, 10*maxy);
	translatePoints(tv);
}

void CSGTerminalNode::initializePoints()
{
// TODO This will probably be a set of circles or elipses
	for(Uint k = 0; k < boundingPoly.size(); ++k)
	{
		Point_2_base tmp(boundingPoly[k].x().to_double(), boundingPoly[k].y().to_double());
		geneList.push_back(W_Point_2(tmp, 0.));
	}
	setAlphaAndWeightsConvex();
	Real sf = HelperNS::RandomGen::instance().randRealInRange(0.6, 0.9);
	scalePoints(sf);
}

void CSGTerminalNode::rotatePoints(const Real theta)
{
	// rotates points about the center of mass of the hull
	Point_2 cg = computeHullCenterOMass();
	
	Real sintheta = sin(theta), costheta = cos(theta);
	std::list<W_Point_2>::iterator plist;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
	{
		Point_2_base tmp = *plist;
		Real x = costheta*(tmp.x() - cg.x().to_double());
		x += sintheta*(tmp.y() - cg.y().to_double());
		x += cg.x().to_double();
		Real y = costheta*(tmp.y() - cg.y().to_double());
		y += sintheta*(tmp.x() - cg.x().to_double());
		y += cg.y().to_double();
		*plist = W_Point_2(Point_2_base(x, y), plist->weight());
	}
}

void CSGTerminalNode::translatePoints(const Vector_2& tp)
{
	Mesh_Vector_2 mtp(tp.x().to_double(), tp.y().to_double());
	std::list<W_Point_2>::iterator plist;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
	{
		Point_2_base tmp = plist->point() + mtp;
		*plist = W_Point_2(tmp, plist->weight());
	}
}

void CSGTerminalNode::scalePoints(const Real sf)
{
	// scales points about the center of mass of the hull
	Point_2 cg = computeHullCenterOMass();
	std::list<W_Point_2>::iterator plist;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
	{
		Real x = sf*plist->x() + (1. - sf)*cg.x().to_double();
		Real y = sf*plist->y() + (1. - sf)*cg.y().to_double();
		*plist = W_Point_2(Point_2_base(x, y), plist->weight());
	}
}

void CSGTerminalNode::generateHull(std::vector<Point_2>& hullVec) const
{
	std::vector<Mesh_Segment_2> segVec;
	decode(segVec);
	std::vector<std::vector<Mesh_Segment_2> > orderedArray;
	GeometryTranslation::orderMeshSegments(segVec, orderedArray);
	for(Uint k = 0; k < orderedArray.size(); ++k)
	{
		vector<Point_2> ptPoly = getPointVecFromPoly(orderedArray[k]);
		hullVec.insert(hullVec.begin(), ptPoly.begin(), ptPoly.end());
	}
}

Point_2 CSGTerminalNode::computeHullCenterOMass()
{
	vector<Point_2> hullVec;
	generateHull(hullVec);

	Vector_2 sum(0.,0.);
	int n = hullVec.size();
	for(int kh = 0; kh < n; ++kh)
		sum = sum + (hullVec[kh] - CGAL::ORIGIN);

	Point_2 cg = CGAL::ORIGIN;
	if(n > 0)
		cg = CGAL::ORIGIN + sum/n;
	
	return cg;
}

void CSGTerminalNode::debugPrint()
{
	cout << "geneList = [";
	for(std::list<W_Point_2>::const_iterator plist = geneList.begin(); plist != geneList.end(); ++plist)
		cout << plist->x() << " " << plist->y() << " " << plist->weight() << endl;
	cout << "];" << endl;
}

unsigned CSGTerminalNode::countPoints() const
{
	return geneList.size();
}

void CSGTerminalNode::printCSGNodeToFile(std::ofstream& popFile) const
{
	popFile << -1 << endl;
	int n = geneList.size();
	popFile << itsGenericMaterial.getNumParameters() << " ";
	for(int k = 0; k < itsGenericMaterial.getNumParameters(); k++)
		popFile << itsGenericMaterial.getParameter(k) << " ";
	vector<Real> rgb = itsGenericMaterial.getPrintColor();
	popFile << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;
	popFile << priority << " " << alpha << endl;
	popFile << n << endl;
	std::list<W_Point_2>::const_iterator plist;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
		popFile << plist->x() << " " << plist->y() << " " << plist->weight() << endl;
	popFile << endl;
}

void CSGTerminalNode::getMPIDataFormat(std::vector<int>& treeVec, std::vector<W_Point_2>& ptVec, 
				 std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const
{
	treeVec.push_back(-1);
	int n = geneList.size();
	treeVec.push_back(n);
	std::list<W_Point_2>::const_iterator plist;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
		ptVec.push_back(*plist);
	priorityVec.push_back(priority);
	priorityVec.push_back(alpha);
	matVec.push_back(itsGenericMaterial);
}

void CSGTerminalNode::getMPIDataFormat(std::vector<int>& treeVec, std::vector<Point_3>& ptVec,
                 std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const
{
}

bool CSGTerminalNode::checkMeshability()
{
	vector<Mesh_Segment_2> tmpAS;
	decode(tmpAS);
	Real a = polyArea(tmpAS);
	Real p = polyPerimeter(tmpAS);
	Real apr = fabs(a)/(p*boundMag);
	if(apr < 0.011 || a == 0. || p == 0. || tmpAS.size() == 0)
	{
		// the min. dimension of this gene is probably too small
		// to reliably mesh this gene, i.e. it's too thin 
		return false;
	}
	return true;
}

Real CSGTerminalNode::polyArea(vector<Mesh_Segment_2>& inAS) const
{
	Real sa = 0.;
	for(Uint k = 0; k < inAS.size(); k++)
	{
		Point_2_base p1 = inAS[k].source(), p2 = inAS[k].target();
		sa += p1.x()*p2.y() - p1.y()*p2.x();
	}
	return 0.5*sa;
}

Real CSGTerminalNode::polyPerimeter(vector<Mesh_Segment_2>& inAS) const
{
	Real sp = 0.;
	for(Uint k = 0; k < inAS.size(); k++)
	{
		Point_2_base p1 = inAS[k].source(), p2 = inAS[k].target();
		sp += sqrt((p1.x() - p2.x())*(p1.x() - p2.x())
			+ (p1.y() - p2.y())*(p1.y() - p2.y()));
	}
	return sp;
}

Real CSGTerminalNode::getMaxPointRadius() const
{
	std::list<W_Point_2>::const_iterator plist;
	Real maxRad = 0.0;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
		maxRad = MAX(maxRad, (plist->point() - CGAL::ORIGIN).squared_length());
	return sqrt(maxRad);
}

bool CSGTerminalNode::isPtInList(Point_2 inPt, std::vector<Point_2>& ptVec)
{
	bool found = false;
	for(unsigned int k = 0; k < ptVec.size() && !found; k++)
	{
		if(inPt == ptVec[k])
			found = true;
	}
	return found;
}

void CSGTerminalNode::boundsCheck(Point_2& inpt) const
{
	// find area coordinates
	int nvp = boundingPoly.size();
	std::vector<FT> areaCoords;
	for(int k = 0; k < nvp; k++)
	{
		int km1 = (k + nvp - 1)%nvp, kp1 = (k + 1)%nvp;
		Point_2 vkm1 = boundingPoly[km1], vk = boundingPoly[k], vkp1 = boundingPoly[kp1];
		FT a = CGAL::area(vkm1, vk, vkp1);
		for(int l = 0; l < nvp; l++)
		{
			if(l != km1 && l != k)
				a *= CGAL::area(inpt, boundingPoly[l], boundingPoly[(l + 1)%nvp]);
		}
		a = MAX(a, 0);
		areaCoords.push_back(a);
	}
	FT sum = 0;
	for(int kv = 0; kv < nvp; kv++)
		sum += areaCoords[kv];

	Vector_2 newPt(0, 0);
	for(int kv = 0; kv < nvp; kv++)
		newPt = newPt + (areaCoords[kv]/sum)*(boundingPoly[kv] - CGAL::ORIGIN);
	inpt = CGAL::ORIGIN + newPt;

	if(useBoundSnap)
		snapToBoundary(inpt);
}

void CSGTerminalNode::boundsCheck(W_Point_2& inpt) const
{
	Point_2_base tmpPt(inpt.x(), inpt.y());
	boundsCheck(tmpPt);
	//Real w = MIN(MAX(inpt.weight(), 0.), weightmax);
	//Real w = MAX(inpt.weight(), 0.);
	Real w = MIN(MAX(inpt.weight(), -weightmax), weightmax);
	inpt = W_Point_2(tmpPt, w);
}

void CSGTerminalNode::boundsCheck(Point_2_base& inpt) const
{
	Point_2 tmpPt(inpt.x(), inpt.y());
	boundsCheck(tmpPt);
	inpt = Point_2_base(tmpPt.x().to_double(), tmpPt.y().to_double());
}

void CSGTerminalNode::snapToBoundary(Point_2& inpt) const
{
	int nvp = boundingPoly.size();
	FT bs2 = boundSnap*boundSnap;
	for(int k = 0; k < nvp; k++)
	{
		int kp1 = (k + 1)%nvp;
		Segment_2 s(boundingPoly[k], boundingPoly[kp1]);
		if(CGAL::squared_distance(inpt, s) < bs2)
		{
			K::Line_2 l = s.supporting_line();
			inpt = l.projection(inpt);
		}
	}
}

W_Point_2 CSGTerminalNode::genRandomPoint() const
{
	Real minx = boundingPoly[0].x().to_double(), maxx = boundingPoly[0].x().to_double(),
		 	 miny = boundingPoly[0].y().to_double(), maxy = boundingPoly[0].y().to_double();
	for(unsigned int k1 = 0; k1 < boundingPoly.size(); k1++)
	{
		minx = MIN(minx, boundingPoly[k1].x().to_double());
		maxx = MAX(maxx, boundingPoly[k1].x().to_double());
		miny = MIN(miny, boundingPoly[k1].y().to_double());
		maxy = MAX(maxy, boundingPoly[k1].y().to_double());
	}
	Real rx = HelperNS::RandomGen::instance().randRealInRange(minx, maxx), 
			 ry = HelperNS::RandomGen::instance().randRealInRange(miny, maxy),
			 rw = HelperNS::RandomGen::instance().randRealInRange(-weightmax, weightmax);
	W_Point_2 newpt(Point_2_base(rx,ry),rw);
	boundsCheck(newpt);
	return newpt;
}

void CSGTerminalNode::setBoundPoly(std::vector<Point_2> inPoly)
{
	boundingPoly = inPoly;
	// find max dimension
	Real minx = boundingPoly[0].x().to_double(), maxx = boundingPoly[0].x().to_double(),
		 miny = boundingPoly[0].y().to_double(), maxy = boundingPoly[0].y().to_double();
	for(unsigned int k1 = 0; k1 < boundingPoly.size(); k1++)
	{
		minx = MIN(minx, boundingPoly[k1].x().to_double());
		maxx = MAX(maxx, boundingPoly[k1].x().to_double());
		miny = MIN(miny, boundingPoly[k1].y().to_double());
		maxy = MAX(maxy, boundingPoly[k1].y().to_double());
	}
	Real xdim = maxx - minx, ydim = maxy - miny;
	boundMag = (xdim + ydim)/2.;
	alphamax = 0.5*boundMag;
	weightmax = 1.*boundMag;
}

void CSGTerminalNode::copyPointList(std::list<W_Point_2>& copyList)
{
    copyList = geneList;
}

void CSGTerminalNode::setPointList(const std::list<W_Point_2>& copyList)
{
    geneList = copyList;
}

bool CSGTerminalNode::isPointInPoly(Point_2& testPt, vector<Point_2>& inPoly) const
{
    if(CGAL::bounded_side_2(inPoly.begin(), inPoly.end(), testPt, K()) == CGAL::ON_BOUNDED_SIDE)
        return true;
    return false;
}

bool CSGTerminalNode::isPointInPoly(Point_2& testPt, vector<Mesh_Segment_2>& inPoly) const
{
	vector<Mesh_Segment_2> segVec;
	decode(segVec);
	std::vector<std::vector<Mesh_Segment_2> > orderedArray;
	GeometryTranslation::orderMeshSegments(segVec, orderedArray);
  for(Uint k = 0; k < orderedArray.size(); ++k)
  {
    vector<Point_2> ptPoly = getPointVecFromPoly(orderedArray[k]);
		if(isPointInPoly(testPt, ptPoly))
			return true;
  }
	return false;
}

void CSGTerminalNode::setAlphaAndWeightsConvex()
{
	alpha = boundMag*10; // Large alpha so we recover convex shapes
	std::list<W_Point_2>::iterator plist1;
  for(plist1 = geneList.begin(); plist1 != geneList.end(); ++plist1)
		*plist1 = W_Point_2(plist1->point(), 0.);
}

void CSGTerminalNode::setVector(std::vector<double>::const_iterator& curPos, const std::vector<double>::const_iterator& endPos)
{
	list<W_Point_2>::iterator lit;
	for(lit = geneList.begin(); lit != geneList.end(); ++lit)
	{
		if(curPos == endPos)
			return;
		double xp = *curPos;
		++curPos;
		if(curPos == endPos)
			return;
		double yp = *curPos;
		++curPos;
		double wp = lit->weight();
		*lit = W_Point_2(Point_2_base(xp, yp), wp);
	}
}

Real CSGTerminalNode::getCurrentMaxWeight() const
{
	Real curMaxW = 0.;
	list<W_Point_2>::const_iterator cit;
	for(cit = geneList.begin(); cit != geneList.end(); ++cit)
		curMaxW = MAX(curMaxW, cit->weight());
	return curMaxW;
}

std::list<W_Point_2>::iterator CSGTerminalNode::getPointIterator(const W_Point_2& findPt)
{
	list<W_Point_2>::iterator foundit = geneList.end(), lit;
	bool found = false;
	for(lit = geneList.begin(); lit != geneList.end() && !found; ++lit)
	{
		if(*lit == findPt)
		{
			found = true;
			foundit = lit;
		}
	}
	return lit;
}
}
