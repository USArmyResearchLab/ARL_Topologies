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

#include "csgterminalnode3d.h"

namespace Topologies{
VerbLevel CSGTerminalNode3D::verbosity = vlSilent;
Polyhedron_3 CSGTerminalNode3D::boundingPoly;

using namespace std;

CSGTerminalNode3D::CSGTerminalNode3D(int curTreeDepth):
	CSGNode(curTreeDepth)
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
		std::size_t rid = HelperNS::RandomGen::instance().randIntInRange((std::size_t) 0, materialList.size() - 1);
		itsGenericMaterial = materialList[rid];
	}
	priority = HelperNS::RandomGen::instance().randRealInRange(2.0, prioritymax);
}

CSGTerminalNode3D::CSGTerminalNode3D(const CSGTerminalNode3D& copy):
	CSGNode(copy),
	itsGenericMaterial(copy.itsGenericMaterial),
	priority(copy.priority)
{
	for(std::list<Point_3>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
		geneList.push_back(*plist);
}

CSGTerminalNode3D::CSGTerminalNode3D(CSGTerminalNode3D&& copy):
	CSGNode(copy)
{
	swap(copy);
}

CSGTerminalNode3D& CSGTerminalNode3D::operator=(CSGTerminalNode3D rhs)
{
	swap(rhs);
	return *this;
}

void CSGTerminalNode3D::swap(CSGTerminalNode3D& arg2)
{
	geneList.swap(arg2.geneList);
	itsGenericMaterial.swap(arg2.itsGenericMaterial);
	std::swap(priority, arg2.priority);
}

CSGTerminalNode3D::CSGTerminalNode3D(const CSGTerminalNode3D& copy, MirrorType inMT):
  CSGNode(copy),
  itsGenericMaterial(copy.itsGenericMaterial),
  priority(copy.priority)
{
	if(inMT == mtNone)
	{
  	for(std::list<Point_3>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
    	geneList.push_back(*plist);
	}
	else if(inMT == mtX)
	{
		for(std::list<Point_3>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
			geneList.push_back(Point_3(-plist->x(), plist->y(), plist->z()));
	}
	else if(inMT == mtY)
	{
		for(std::list<Point_3>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
      geneList.push_back(Point_3(plist->x(), -plist->y(), plist->z()));
	}
	else if(inMT == mtZ)
	{
		for(std::list<Point_3>::const_iterator plist = (copy.geneList).begin(); plist != (copy.geneList).end(); ++plist)
			geneList.push_back(Point_3(plist->x(), plist->y(), -plist->z()));
	}
}

CSGTerminalNode3D::CSGTerminalNode3D(int curTreeDepth, std::ifstream& saveFile):
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
	saveFile >> numPoints;
	for(int k = 0; k < numPoints; k++)
	{
		Real x, y, z;
		saveFile >> x >> y >> z;
		Point_3 newPoint(x, y, z);
		boundsCheck(newPoint);
		geneList.push_back(newPoint);
	}
}

CSGTerminalNode3D::CSGTerminalNode3D(int curTreeDepth, int numMatParams, int& ptVecPos, int& matPos, int& prPos, int& treeVecPos, 
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
	for(int k = 0; k < numPoints; k++)
	{
		Real x, y, z;
		x = ptArray[ptVecPos++];
		y = ptArray[ptVecPos++];
		z = ptArray[ptVecPos++];
		Point_3 newPoint(x, y, z);
		geneList.push_back(newPoint);
	}
}

CSGTerminalNode3D::CSGTerminalNode3D(int curTreeDepth, int numMatParams, std::vector<double>::const_iterator& ptIt,
                    std::vector<double>::const_iterator& matIt, std::vector<double>::const_iterator& prIt,
                    std::vector<int>::const_iterator& treeIt) :
	CSGNode(curTreeDepth)
{
	std::vector<Real> constParams;
	for(int k = 0; k < numMatParams; k++)
	{
		constParams.push_back(*matIt);
		matIt++;
	}
	itsGenericMaterial = GenericMaterial(constParams);
	priority = *prIt;
	prIt++;
	unsigned numPoints = *treeIt;
	++treeIt;
	for(int k = 0; k < numPoints; k++)
	{
		Real x, y, z;
		x = *ptIt;
		++ptIt;
		y = *ptIt;
		++ptIt;
		z = *ptIt;
		++ptIt;
		Point_3 newPoint(x, y, z);
		geneList.push_back(newPoint);
	}
}

CSGTerminalNode3D::CSGTerminalNode3D(std::list<Point_3>& ptList):
	CSGNode(0)
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
}

CSGTerminalNode3D::~CSGTerminalNode3D()
{
}

bool CSGTerminalNode3D::decodeMaterial(Point_3& testPt, GenericMaterial& theGM) const
{
	Polyhedron_3 tmpCH;
	generateHull(tmpCH);
	if(isPointInPoly(testPt, tmpCH))
	{
		theGM = itsGenericMaterial;
		return true;
	}
	return false;
}

bool CSGTerminalNode3D::decodeMaterial(Point_3& testPt, int& theGM) const
{
	theGM++;
    Polyhedron_3 tmpCH;
    generateHull(tmpCH);
    if(isPointInPoly(testPt, tmpCH))
        return true;
    return false;
}

void CSGTerminalNode3D::setMaterialList(std::vector<GenericMaterial>& theGMVec) const
{
	theGMVec.push_back(itsGenericMaterial);
}

void CSGTerminalNode3D::decode(Nef_polyhedron_3& inNP) const
{
	Polyhedron_3 newP;
	generateHull(newP);
    if(newP.is_closed())
	{
		// exclude boundaries
		inNP = Nef_polyhedron_3(newP);
		inNP = inNP.interior();
	}
}

void CSGTerminalNode3D::decodeMask(Nef_polyhedron_2& inNP) const
{
}

bool CSGTerminalNode3D::decodeMaterial(Point_2& testPt, GenericMaterial& theGM) const
{
    return false;
}

void CSGTerminalNode3D::decode(Nef_polyhedron_2& inNP) const
{
}

void CSGTerminalNode3D::decodeMask(Nef_polyhedron_3& inNP) const
{
    Polyhedron_3 newP;
    generateHull(newP);
    if(newP.is_closed())
        inNP = Nef_polyhedron_3(newP);
}

void CSGTerminalNode3D::movePointsOB()
{
	Vertex_iterator_3 vit = boundingPoly.vertices_begin();
	Real maxx = vit->point().x().to_double(), maxy = vit->point().y().to_double(),
		maxz = vit->point().z().to_double();
	for(++vit; vit != boundingPoly.vertices_end(); ++vit)
	{
        maxx = MAX(maxx, vit->point().x().to_double());
        maxy = MAX(maxy, vit->point().y().to_double());
        maxz = MAX(maxz, vit->point().z().to_double());
    }
	Vector_3 tv(10*maxx, 10*maxy, 10*maxz);
	Real sf = 0.001*maxx;
	scalePoints(sf);
	translatePoints(tv);
}

void CSGTerminalNode3D::getMPIDataFormat(std::vector<int>& treeVec, std::vector<W_Point_2>& ptVec,
                  std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const
{
}

void CSGTerminalNode3D::rotatePoints(const Real alpha, const Real beta, const Real gamma)
{
    // rotates points about the center of mass of the hull
    Vector_3 cg = computeHullCenterOMass() - CGAL::ORIGIN;

    Real sinvar = sin(alpha), cosvar = cos(alpha);
    list<Point_3>::iterator plist;
    for(plist = geneList.begin(); plist != geneList.end(); ++plist)
    {
        *plist = *plist - cg;
        Real x = (*plist)[0].to_double(), y = (*plist)[1].to_double(), z = (*plist)[2].to_double();
        Real yn = cosvar*y + sinvar*z;
        Real zn = cosvar*z - sinvar*y;
        *plist = Point_3(x, yn, zn);
        *plist = *plist + cg;
    }

    sinvar = sin(beta);
    cosvar = cos(beta);
    for(plist = geneList.begin(); plist != geneList.end(); ++plist)
    {
        *plist = *plist - cg;
        Real x = (*plist)[0].to_double(), y = (*plist)[1].to_double(), z = (*plist)[2].to_double();
        Real xn = cosvar*x - sinvar*z;
        Real zn = cosvar*z + sinvar*x;
        *plist = Point_3(xn, y, zn);
        *plist = *plist + cg;
    }

    sinvar = sin(gamma);
    cosvar = cos(gamma);
    for(plist = geneList.begin(); plist != geneList.end(); ++plist)
    {
        *plist = *plist - cg;
        Real x = (*plist)[0].to_double(), y = (*plist)[1].to_double(), z = (*plist)[2].to_double();
        Real xn = cosvar*x + sinvar*y;
        Real yn = cosvar*y - sinvar*x;
        *plist = Point_3(xn, yn, z);
        *plist = *plist + cg;
    }
}

void CSGTerminalNode3D::translatePoints(const Vector_3& tp)
{
    list<Point_3>::iterator plist;
    for(plist = geneList.begin(); plist != geneList.end(); ++plist)
        *plist = *plist + tp;
}

void CSGTerminalNode3D::scalePoints(const Real sf)
{
    // scales points about the center of mass of the hull
    Point_3 cg = computeHullCenterOMass();

    list<Point_3>::iterator plist;
    for(plist = geneList.begin(); plist != geneList.end(); ++plist)
    {
        Real x, y, z;
        x = sf*(*plist)[0].to_double() + (1 - sf)*cg[0].to_double();
        y = sf*(*plist)[1].to_double() + (1 - sf)*cg[1].to_double();
        z = sf*(*plist)[2].to_double() + (1 - sf)*cg[2].to_double();
        *plist = Point_3(x, y, z);
    }
}

void CSGTerminalNode3D::generateHull(Polyhedron_3& P) const
{
	// compute convex hull
	std::list<Point_3> ptVec;
	std::list<Point_3>::const_iterator plist;

	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
	{
		// check to make sure this point is not already in list;
		// coincidental points crash the mesher apparently
		bool found = false;
		std::list<Point_3>::const_iterator plist2;
		for(plist2 = ptVec.begin(); plist2 != ptVec.end() && !found; ++plist2)
		{
			if(sqrt((*plist - *plist2).squared_length().to_double()) < boundMag*1.e-8)
				found = true;
		}
		if(!found)
			ptVec.push_back(*plist);
	}
	CGAL::Object ch_object;
	CGAL::convex_hull_3(ptVec.begin(), ptVec.end(), ch_object);
	CGAL::assign(P, ch_object);
}

void CSGTerminalNode3D::generateHull(vector<Point_3>& hullVec) const
{
	// compute convex hull
	std::list<Point_3> ptVec;
	std::list<Point_3>::const_iterator plist;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
		ptVec.push_back(*plist);
	Polyhedron_3 P;
	CGAL::Object ch_object;
	CGAL::convex_hull_3(geneList.begin(), geneList.end(), ch_object);
	CGAL::assign(P, ch_object);
	hullVec.clear();
	for(Vertex_iterator_3 v = P.vertices_begin(); v != P.vertices_end(); ++v)
		hullVec.push_back(v->point());
}

Point_3 CSGTerminalNode3D::computeHullCenterOMass()
{
	if(verbosity == vlDebug) cout << "getting center of mass" << endl;
	std::vector<Point_3> hullVec;
	generateHull(hullVec);
	
	Vector_3 sum(0.,0., 0.);
	int n = hullVec.size();
	for(int kh = 0; kh < n; ++kh)
		sum = sum + (hullVec[kh] - CGAL::ORIGIN);
	Point_3 cg = CGAL::ORIGIN + sum/static_cast<Real>(n);
	if(verbosity == vlDebug) cout << "done c.o.m." << endl;
	return cg;
}

void CSGTerminalNode3D::debugPrint()
{
	cout << "geneList = [";
	for(std::list<Point_3>::const_iterator plist = geneList.begin(); plist != geneList.end(); ++plist)
	{
		cout << (*plist).x().to_double() << " " << (*plist).y().to_double() << " " << (*plist).z().to_double()<< endl;
	}
	cout << "];" << endl;
}

unsigned CSGTerminalNode3D::countPoints() const
{
	return geneList.size();
}

void CSGTerminalNode3D::printCSGNodeToFile(std::ofstream& popFile) const
{
	popFile << -1 << endl;
	int n = geneList.size();
	popFile << itsGenericMaterial.getNumParameters() << " ";
	for(int k = 0; k < itsGenericMaterial.getNumParameters(); k++)
		popFile << itsGenericMaterial.getParameter(k) << " ";
	vector<Real> rgb = itsGenericMaterial.getPrintColor();
	popFile << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;
	popFile << priority << endl;
	popFile << n << endl;
	std::list<Point_3>::const_iterator plist;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
		popFile << (*plist).x().to_double() << " " << (*plist).y().to_double() << " " << (*plist).z().to_double() << endl;
	popFile << endl;
}

void CSGTerminalNode3D::getMPIDataFormat(std::vector<int>& treeVec, std::vector<Point_3>& ptVec, 
				 std::vector<GenericMaterial>& matVec, std::vector<Real>& priorityVec) const
{
	treeVec.push_back(-1);
	int n = geneList.size();
	treeVec.push_back(n);
	std::list<Point_3>::const_iterator plist;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
		ptVec.push_back(*plist);
	priorityVec.push_back(priority);
	matVec.push_back(itsGenericMaterial);
}

bool CSGTerminalNode3D::checkMeshability()
{
	return true;
}

Real CSGTerminalNode3D::getMaxPointRadius() const
{
	std::list<Point_3>::const_iterator plist;
	Real maxRad = 0.0;
	for(plist = geneList.begin(); plist != geneList.end(); ++plist)
		maxRad = MAX(maxRad, ((*plist) - CGAL::ORIGIN).squared_length().to_double());
	return sqrt(maxRad);
}

void CSGTerminalNode3D::boundsCheck(Point_3& inpt)
{
	if(!isPointInPoly(inpt, boundingPoly))
	{
		if(verbosity == vlDebug) cout << "moving point to poly boundary" << endl;
		// point is out of bounds, project onto a defining plane
		Facet_iterator_3 fit = boundingPoly.facets_begin();
		Facet_iterator_3 closestFit = fit;
		Plane_3 closestPlane = getFacetPlane(fit);
		FT pdist = CGAL::squared_distance(closestPlane, inpt);
		for(fit++; fit != boundingPoly.facets_end(); fit++)
		{
			Plane_3 tmpp = getFacetPlane(fit);
			FT tmpdist = CGAL::squared_distance(tmpp, inpt);
			if(tmpdist < pdist)
			{
				pdist = tmpdist;
				closestFit = fit;
			}
		}
		inpt = getFacetPlane(closestFit).projection(inpt);
		if(verbosity == vlDebug) cout << "done" << endl;
	}
}

Plane_3 CSGTerminalNode3D::getFacetPlane(Facet_iterator_3& fit) const
{
	std::vector<Point_3> planeVec;
    HE_circulator_3 h = fit->facet_begin();
    do
    {
        planeVec.push_back(h->vertex()->point());
        ++h;
    } while(h != fit->facet_begin());
    Plane_3 curplane(planeVec[0], planeVec[1], planeVec[2]);
	return curplane;
}

Point_3 CSGTerminalNode3D::genRandomPoint()
{
	std::vector<Real> lambda;
	Real sum = 0.;
	for(Uint kv = 0; kv < boundingPoly.size_of_vertices(); kv++)
	{
		lambda.push_back(HelperNS::RandomGen::instance().randRealInRange<double>());
		sum += lambda[kv];
	}
	Vector_3 newPt(0., 0., 0.);
	Point_iterator_3 pit;
	std::vector<Real>::const_iterator lit = lambda.begin();
	int kc = 0;
	for(pit = boundingPoly.points_begin(); pit != boundingPoly.points_end(); ++pit)
	{
		newPt = newPt + ((*lit)/sum)*(*pit - CGAL::ORIGIN);
		lit++;
	}
	return CGAL::ORIGIN + newPt;
}

Point_3 CSGTerminalNode3D::genRandomPointInHull()
{
	vector<Point_3> hullVec;
	generateHull(hullVec);
    std::vector<Real> lambda;
    Real sum = 0.;
    int nvp = hullVec.size();
    for(int kv = 0; kv < nvp; kv++)
    {
        lambda.push_back(HelperNS::RandomGen::instance().randRealInRange<double>());
        sum += lambda[kv];
    }
    Vector_3 newPt(0., 0., 0.);
    for(int kv = 0; kv < nvp; kv++)
        newPt = newPt + (lambda[kv]/sum)*(hullVec[kv] - CGAL::ORIGIN);
    return CGAL::ORIGIN + newPt;
}

void CSGTerminalNode3D::setBoundPoly(Polyhedron_3& inPoly)
{
	boundingPoly = inPoly;
	if(verbosity == vlDebug) cout << "bound poly set: " << endl;
	cout << boundingPoly << endl;
	// find max dimension
	Vertex_iterator_3 vit = boundingPoly.vertices_begin();
	Real minx = vit->point().x().to_double(), maxx = vit->point().x().to_double(),
		 miny = vit->point().y().to_double(), maxy = vit->point().y().to_double(),
		 minz = vit->point().z().to_double(), maxz = vit->point().z().to_double();
	for(++vit; vit != boundingPoly.vertices_end(); ++vit)
	{
		minx = MIN(minx, vit->point().x().to_double());
		maxx = MAX(maxx, vit->point().x().to_double());
		miny = MIN(miny, vit->point().y().to_double());
		maxy = MAX(maxy, vit->point().y().to_double());
		minz = MIN(minz, vit->point().z().to_double());
        maxz = MAX(maxz, vit->point().z().to_double());
	}
 
	Real xdim = maxx - minx, ydim = maxy - miny, zdim = maxz - minz;
	boundMag = (xdim + ydim + zdim)/3.;
}

CGAL::Bbox_3 CSGTerminalNode3D::getBoundBox()
{
    // find max dimension
    Vertex_iterator_3 vit = boundingPoly.vertices_begin();
    Real minx = vit->point().x().to_double(), maxx = vit->point().x().to_double(),
         miny = vit->point().y().to_double(), maxy = vit->point().y().to_double(),
         minz = vit->point().z().to_double(), maxz = vit->point().z().to_double();
    for(++vit; vit != boundingPoly.vertices_end(); ++vit)
    {
        minx = MIN(minx, vit->point().x().to_double());
        maxx = MAX(maxx, vit->point().x().to_double());
        miny = MIN(miny, vit->point().y().to_double());
        maxy = MAX(maxy, vit->point().y().to_double());
        minz = MIN(minz, vit->point().z().to_double());
        maxz = MAX(maxz, vit->point().z().to_double());
    }
	Real ff = 0.01;
	Real xf = ff*(maxx - minx), yf = ff*(maxy - miny), zf = ff*(maxz - minz);
    return CGAL::Bbox_3(minx - xf, miny - yf, minz - zf, maxx + xf, maxy + yf, maxz + zf);
}

bool CSGTerminalNode3D::isPointInBoundary(const Point_3& testPt)
{
	return isPointInPoly(testPt, boundingPoly);
}


void CSGTerminalNode3D::copyPointList(std::list<Point_3>& copyList)
{
    copyList = geneList;
}

void CSGTerminalNode3D::setPointList(const std::list<Point_3>& copyList)
{
    geneList = copyList;
}

bool CSGTerminalNode3D::isPointInPoly(const Point_3& tstPt, Polyhedron_3& poly)
{
	Facet_iterator_3 fit;
	bool first = true, possign, allsame = true;
	for(fit = poly.facets_begin(); fit != poly.facets_end(); fit++)
	{
		HE_circulator_3 hec = fit->facet_begin();
		vector<Point_3> faceVec;
		for(Uint k = 0; k < 3; k++)
		{
			faceVec.push_back(hec->vertex()->point());
			++hec;
		}
		bool sgn = CGAL::volume(faceVec[0], faceVec[1], faceVec[2], tstPt) >= 0;
		if(first)
		{
			possign = sgn;
			first = false;
		}
		else 
			allsame &= (possign == sgn);
	}
	return allsame;
}

void CSGTerminalNode3D::setVector(std::vector<double>::const_iterator& curPos, const std::vector<double>::const_iterator& endPos)
{
	list<Point_3>::iterator lit;
	Uint k = 1;
	for(lit = geneList.begin(); lit != geneList.end(); ++lit)
	{
		Real xp, yp, zp;
		if(curPos == endPos);
			return;
		xp = *curPos;
		++curPos;
		if(curPos == endPos);
			return;
		yp = *curPos;
		++curPos;
		if(curPos == endPos);
			return;
		zp = *curPos;
		++curPos;
		*lit = Point_3(xp, yp, zp);
	}
}

void CSGTerminalNode3D::getLocalOptFormat(std::list<W_Point_2>& ptList) const
{
}
}

