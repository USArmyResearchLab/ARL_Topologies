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

#include "csgfunctionnode.h"
#include "csgterminalnode.h"
#include "csgterminalnode3d.h"
#include <memory>

namespace Topologies{
using std::vector;
using std::cout;
using std::endl;

CSGFunctionNode::CSGFunctionNode(int curTreeDepth, unsigned numNodesX, unsigned numNodesY, unsigned numPts, 
																 double width, double height):
	CSGNode(curTreeDepth)
{	
	// Function type defaults to union for now
/*
	if(depthInTree == 1)
	{
		theFT = ftSubtraction;
		children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(curTreeDepth + 1, 16, Point_2_base(width*0.5, height*0.5), 0.5, 0.5)));
//		children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(depthInTree, numPts, Point_2_base(width*0.5, height*0.5), 
//																																		width*2, height*2)));
		children.push_back(std::unique_ptr<CSGNode>(new CSGFunctionNode(depthInTree, numNodesX, numNodesY, numPts, width, height)));
	}
	else
	{*/
		theFT = ftUnion;
		if(numNodesX*numNodesY > 0)
			createChildren(numNodesX, numNodesY, numPts, width, height);
//	}
}

CSGFunctionNode::CSGFunctionNode(int curTreeDepth, const std::vector<std::list<W_Point_2>>& ptListVec) :
	CSGNode(curTreeDepth)
{
	theFT = ftUnion;
	assert(theDim == dt2d);
	for(auto it = ptListVec.begin(); it != ptListVec.end(); ++it)
		children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(depthInTree, *it)));
}

CSGFunctionNode::CSGFunctionNode(const CSGFunctionNode& copy):
	CSGNode(copy),
	theFT(copy.theFT)
{
	int nc = copy.children.size();
	for(int kc = 0; kc < nc; kc++)
	{
		CSGNode* ngChild = copy.children[kc].get();
		if(ngChild->getNumChildren() == 0)
		{
			if(theDim == dt2d)
				children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(*dynamic_cast<CSGTerminalNode*>(ngChild))));
			else
				children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode3D(*dynamic_cast<CSGTerminalNode3D*>(ngChild))));
		}
		else
			children.push_back(std::unique_ptr<CSGNode>(new CSGFunctionNode(*dynamic_cast<CSGFunctionNode*>(ngChild))));
	}
}

CSGFunctionNode::CSGFunctionNode(CSGFunctionNode && copy) :
	CSGNode(copy)
{
	swap(copy);
}

CSGFunctionNode& CSGFunctionNode::operator=(CSGFunctionNode copy)
{
	swap(copy);
	return *this;
}

void CSGFunctionNode::swap(CSGFunctionNode& arg2)
{
	CSGNode::swap(arg2);
	std::swap(theFT, arg2.theFT);
	children.swap(arg2.children);
}

CSGFunctionNode::CSGFunctionNode(const CSGFunctionNode& copy, MirrorType inMT):
  theFT(copy.theFT),
  CSGNode(copy)
{
  int nc = copy.children.size();
  for(int kc = 0; kc < nc; kc++)
  {
    CSGNode* ngChild = copy.children[kc].get();
    if(ngChild->getNumChildren() == 0)
    {
      if(theDim == dt2d)
        children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(*dynamic_cast<CSGTerminalNode*>(ngChild), inMT)));
      else
        children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode3D(*dynamic_cast<CSGTerminalNode3D*>(ngChild), inMT)));
    }
    else
      children.push_back(std::unique_ptr<CSGNode>(new CSGFunctionNode(*dynamic_cast<CSGFunctionNode*>(ngChild), inMT)));
  }
}

CSGFunctionNode::CSGFunctionNode(int curTreeDepth, int functionType, std::ifstream& saveFile):
	theFT(static_cast<FunctionType>(functionType)),
	CSGNode(curTreeDepth)
{
	int numKids;
	saveFile >> numKids;
	for(int k = 0; k < numKids; k++)
	{
		int nodeType;
		saveFile >> nodeType;
		if(nodeType >= 0)
			children.push_back(std::unique_ptr<CSGNode>(new CSGFunctionNode(depthInTree, nodeType, saveFile)));
		else
		{
			if(theDim == dt2d)
				children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(depthInTree, saveFile)));
			else
				children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode3D(depthInTree, saveFile)));
		}
	}
}

CSGFunctionNode::CSGFunctionNode(int curTreeDepth, int numMatParams, int functionType, int& ptVecPos, int& matPos, int& prPos,
							int& treeVecPos, Real* ptArray, Real* matArray, Real* prArray, int* treeArray):
	theFT(static_cast<FunctionType>(functionType)),
	CSGNode(curTreeDepth)
{
	int numKids = treeArray[treeVecPos++];
	for(int k = 0; k < numKids; k++)
	{
		int nodeType = treeArray[treeVecPos++];
		if(nodeType >= 0)
			children.push_back(std::unique_ptr<CSGNode>(new CSGFunctionNode(depthInTree, numMatParams, nodeType, ptVecPos, 
																																			matPos, prPos, treeVecPos, ptArray, matArray, 
																																			prArray, treeArray)));
		else
		{
			if(theDim == dt2d)
				children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(depthInTree, numMatParams, ptVecPos, matPos, prPos,
																																			 treeVecPos, ptArray, matArray, prArray, treeArray)));
			else	
				children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode3D(depthInTree, numMatParams, ptVecPos, matPos, 
																																					prPos, treeVecPos, ptArray, matArray, prArray, 
																																					treeArray)));
		}
	}
}

CSGFunctionNode::CSGFunctionNode(int curTreeDepth, int functionType, std::vector<int>::const_iterator& nmIt, 
																	std::vector<double>::const_iterator& ptIt, std::vector<double>::const_iterator& matIt, 
																	std::vector<double>::const_iterator& prIt, std::vector<int>::const_iterator& treeIt) :
	theFT(static_cast<FunctionType>(functionType)),
	CSGNode(curTreeDepth)
{
	int numKids = *treeIt;
	++treeIt;
	for(int k = 0; k < numKids; k++)
	{
		int nodeType = *treeIt;
		++treeIt;
		if(nodeType >= 0)
			children.push_back(std::unique_ptr<CSGNode>(new CSGFunctionNode(depthInTree, nodeType, nmIt, ptIt, matIt, prIt, treeIt)));
		else
		{
			int numMatParams = *nmIt;
			++nmIt;
			if(theDim == dt2d)
				children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(depthInTree, numMatParams, ptIt, matIt, 
																																								prIt, treeIt)));
			else
				children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode3D(depthInTree, numMatParams, ptIt, matIt, 
																																									prIt, treeIt)));
		}
	}
}

CSGFunctionNode::CSGFunctionNode(int curTreeDepth, FunctionType inFT, CSGNode* inParent, 
																std::vector<std::unique_ptr<CSGNode> >& inKids):
	CSGNode(curTreeDepth),
	theFT(inFT)
{
	parent = inParent;
	children = std::move(inKids);
}

CSGFunctionNode::~CSGFunctionNode()
{
}

void CSGFunctionNode::createChildren(unsigned numNodesX, unsigned numNodesY, unsigned numPts, double width, double height)
{
	for(unsigned ky = 0; ky < numNodesY; ++ky)
	{
		double dy = height/((double)numNodesY + 1.);
		unsigned curNumX = numNodesX;
		if((ky % 2) == 1 && numNodesX > 0 && ky != (numNodesY - 1))	
			--curNumX;
		for(unsigned kx = 0; kx < curNumX; ++kx)
		{
			double dx = width/((double)curNumX + 1.);
			Point_2_base center(((double)kx + 1.)*dx, ((double)ky + 1.)*dy);
// TODO:: Add 3d
			children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(depthInTree, numPts, center, 0.3*dx, 0.3*dy)));
		}
	}
}

void CSGFunctionNode::createChildren(unsigned numNodes)
{
	int nc = 2;
	numNodes += nc;
	for(int kc = 0; kc < nc; kc++)
	{
		if(theDim == dt2d)
			children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode(depthInTree)));
		else
			children.push_back(std::unique_ptr<CSGNode>(new CSGTerminalNode3D(depthInTree)));
	}
}

bool CSGFunctionNode::decodeMaterial(Point_2& testPt, GenericMaterial& theGM) const
{
	vector<bool> bvec;
	vector<GenericMaterial> gmvec;
	bool allTrue = true, oneTrue = false;
	int nc = children.size();

	for(int kc = 0; kc < nc; kc++)
	{
		GenericMaterial tmp = theGM;
		bvec.push_back(children[kc]->decodeMaterial(testPt, tmp));
		allTrue &= bvec[kc];
		oneTrue |= bvec[kc];
		gmvec.push_back(tmp);
	}
	if(theFT == ftUnion)
	{
		if(allTrue)
			theGM = gmvec[nc-1];
		else if(oneTrue)
		{
			for(int kc = 0; kc < nc; kc++)
			{
				if(bvec[kc])
					theGM = gmvec[kc];
			}
		}
		return allTrue || oneTrue;
	}
	else if(theFT == ftSubtraction)
	{
		bool inNull = true;
		for(int kc = 1; kc < nc; kc++)
			inNull &= bvec[kc];
		if(bvec[0] && !inNull)
		{
			theGM = gmvec[0];
			return true;
		}
	}
	return false;
}

void CSGFunctionNode::decode(Nef_polyhedron_2& inNP) const
{
	vector<Nef_polyhedron_2> sdvec_exc;
	Nef_polyhedron_2 rop;
	int nc = children.size();
	for(int kc = 0; kc < nc; kc++)
	{
		Nef_polyhedron_2 tmpNP1, tmpNP2;
		if(kc > 0)
		{
			children[kc]->decodeMask(tmpNP1);
			rop += tmpNP1;
		}
		children[kc]->decode(tmpNP2);
		sdvec_exc.push_back(tmpNP2);
	}
	if(theFT == ftUnion)
	{
		if(theDT == dtBoolean)
			unionFunction(sdvec_exc, inNP);
		else
			overlapFunction(sdvec_exc, rop, inNP);
	}
	else if(theFT == ftSubtraction)
		subtractionFunction(sdvec_exc, rop, inNP);
}

void CSGFunctionNode::setMaterialList(std::vector<GenericMaterial>& theGMVec) const
{
	Uint nc = children.size();
	for(Uint kc = 0; kc < nc; kc++)
    	children[kc]->setMaterialList(theGMVec);
}

void CSGFunctionNode::decodeMask(Nef_polyhedron_2& inNP) const
{
	// returns the edge-included mask of the desired inhomogeneous geometry
    vector<Nef_polyhedron_2> sdvec;
    int nc = children.size();
    for(int kc = 0; kc < nc; kc++)
    {
        Nef_polyhedron_2 tmpNP;
        children[kc]->decodeMask(tmpNP);
		if(kc > 0 && theFT == ftSubtraction)
			children[kc]->decode(tmpNP);
        sdvec.push_back(tmpNP);
    }
    if(theFT == ftUnion)
        unionFunction(sdvec, inNP);
    else if(theFT == ftSubtraction)
        subtractionFunction(sdvec, inNP);
}

bool CSGFunctionNode::decodeMaterial(Point_3& testPt, GenericMaterial& theGM) const
{
    vector<bool> bvec;
    vector<GenericMaterial> gmvec;
    bool allTrue = true, oneTrue = false;
    int nc = children.size();

    for(int kc = 0; kc < nc; kc++)
    {
        GenericMaterial tmp = theGM;
        bvec.push_back(children[kc]->decodeMaterial(testPt, tmp));
        allTrue &= bvec[kc];
        oneTrue |= bvec[kc];
        gmvec.push_back(tmp);
    }
    if(theFT == ftUnion)
    {
        if(allTrue)
            theGM = gmvec[nc-1];
        else if(oneTrue)
        {
            for(int kc = 0; kc < nc; kc++)
            {
                if(bvec[kc])
                    theGM = gmvec[kc];
            }
        }
        return allTrue || oneTrue;
    }
    else if(theFT == ftSubtraction)
    {
        bool inNull = true;
        for(int kc = 1; kc < nc; kc++)
            inNull &= bvec[kc];
        if(bvec[0] && !inNull)
        {
            theGM = gmvec[0];
            return true;
        }
    }
    return false;
}

bool CSGFunctionNode::decodeMaterial(Point_3& testPt, int& theGM) const
{
    vector<bool> bvec;
    vector<int> gmvec;
    bool allTrue = true, oneTrue = false;
    int nc = children.size();

    for(int kc = 0; kc < nc; kc++)
    {
        bvec.push_back(children[kc]->decodeMaterial(testPt, theGM));
        allTrue &= bvec[kc];
        oneTrue |= bvec[kc];
        gmvec.push_back(theGM);
    }
    if(theFT == ftUnion)
    {
        if(allTrue)
            theGM = gmvec[nc-1];
        else if(oneTrue)
        {
            for(int kc = 0; kc < nc; kc++)
            {
                if(bvec[kc])
                    theGM = gmvec[kc];
            }
        }
        return allTrue || oneTrue;
    }
    else if(theFT == ftSubtraction)
    {
        bool inNull = true;
        for(int kc = 1; kc < nc; kc++)
            inNull &= bvec[kc];
        if(bvec[0] && !inNull)
        {
            theGM = gmvec[0];
            return true;
        }
    }
    return false;
}

void CSGFunctionNode::decode(Nef_polyhedron_3& inNP) const
{
    vector<Nef_polyhedron_3> sdvec_exc;
    Nef_polyhedron_3 rop;
    int nc = children.size();
    for(int kc = 0; kc < nc; kc++)
    {
        Nef_polyhedron_3 tmpNP1, tmpNP2;
        if(kc > 0)
        {
            children[kc]->decodeMask(tmpNP1);
            rop += tmpNP1;
        }
        children[kc]->decode(tmpNP2);
        sdvec_exc.push_back(tmpNP2);
    }
    if(theFT == ftUnion)
    {
        if(theDT == dtBoolean)
            unionFunction(sdvec_exc, inNP);
        else
            overlapFunction(sdvec_exc, rop, inNP);
    }
    else if(theFT == ftSubtraction)
        subtractionFunction(sdvec_exc, rop, inNP);
}

void CSGFunctionNode::decodeMask(Nef_polyhedron_3& inNP) const
{
    // returns the edge-included mask of the desired inhomogeneous geometry
    vector<Nef_polyhedron_3> sdvec;
    int nc = children.size();
    for(int kc = 0; kc < nc; kc++)
    {
        Nef_polyhedron_3 tmpNP;
        children[kc]->decodeMask(tmpNP);
        if(kc > 0 && theFT == ftSubtraction)
            children[kc]->decode(tmpNP);
        sdvec.push_back(tmpNP);
    }
    if(theFT == ftUnion)
        unionFunction(sdvec, inNP);
    else if(theFT == ftSubtraction)
        subtractionFunction(sdvec, inNP);
}

void CSGFunctionNode::movePointsOB()
{
    int nc = children.size();
    for(int kc = 0; kc < nc; kc++)
		children[kc]->movePointsOB();
}

void CSGFunctionNode::unionFunction(vector<Nef_polyhedron_2>& sdvec, Nef_polyhedron_2& inNP) const
{
    int ns = sdvec.size();
    inNP = sdvec[0];
    for(int ks = 1; ks < ns; ks++)
        inNP += sdvec[ks];
}

void CSGFunctionNode::subtractionFunction(vector<Nef_polyhedron_2>& sdvec, Nef_polyhedron_2& inNP) const
{
    int ns = sdvec.size();
    inNP = sdvec[0];
    for(int ks = 1; ks < ns; ks++)
        inNP -= sdvec[ks];
}

void CSGFunctionNode::subtractionFunction(vector<Nef_polyhedron_2>& sdvec, Nef_polyhedron_2& rop, Nef_polyhedron_2& inNP) const
{
    int ns = sdvec.size();
    inNP = sdvec[0];
    inNP -= rop;
}

void CSGFunctionNode::overlapFunction(vector<Nef_polyhedron_2>& sdvec_exc, Nef_polyhedron_2& rop, Nef_polyhedron_2& inNP) const
{
	int ns = sdvec_exc.size();
	inNP = sdvec_exc[0];
	for(int ks = 1; ks < ns; ks++)
	{
		Nef_polyhedron_2 tmp = inNP - rop; 
		inNP = tmp + sdvec_exc[ks];
	}
}

void CSGFunctionNode::unionFunction(vector<Nef_polyhedron_3>& sdvec, Nef_polyhedron_3& inNP) const
{
    int ns = sdvec.size();
    inNP = sdvec[0];
    for(int ks = 1; ks < ns; ks++)
        inNP += sdvec[ks];
}

void CSGFunctionNode::subtractionFunction(vector<Nef_polyhedron_3>& sdvec, Nef_polyhedron_3& inNP) const
{
    int ns = sdvec.size();
    inNP = sdvec[0];
    for(int ks = 1; ks < ns; ks++)
        inNP -= sdvec[ks];
}

void CSGFunctionNode::subtractionFunction(vector<Nef_polyhedron_3>& sdvec, Nef_polyhedron_3& rop, Nef_polyhedron_3& inNP) const
{
    int ns = sdvec.size();
    inNP = sdvec[0];
    inNP -= rop;
}

void CSGFunctionNode::overlapFunction(vector<Nef_polyhedron_3>& sdvec_exc, Nef_polyhedron_3& rop, Nef_polyhedron_3& inNP) const
{
    int ns = sdvec_exc.size();
    inNP = sdvec_exc[0];
    for(int ks = 1; ks < ns; ks++)
    {
        Nef_polyhedron_3 tmp = inNP - rop;
        inNP = tmp + sdvec_exc[ks];
    }
}

unsigned CSGFunctionNode::countPoints() const
{
	unsigned numPointsInSubTree = 0;
	unsigned nc = children.size();
	for(int kc = 0; kc < nc; kc++)
		numPointsInSubTree += children[kc]->countPoints();
	return numPointsInSubTree;
}

unsigned CSGFunctionNode::countNodes() const
{
	unsigned numNodesInSubTree = 0;
	unsigned nc = children.size();
	for(int kc = 0; kc < nc; kc++)
		numNodesInSubTree += children[kc]->countPoints();
	return numNodesInSubTree;
}

void CSGFunctionNode::printCSGNodeToFile(std::ofstream& popFile) const
{
	popFile << static_cast<int>(theFT) << endl;
	int n = children.size();
	popFile << n << endl;
	for(int k = 0; k < n; k++)
		children[k]->printCSGNodeToFile(popFile);
}

void CSGFunctionNode::getMPIDataFormat(vector<int>& treeVec, vector<W_Point_2>& ptVec, 
					vector<GenericMaterial>& matVec, vector<Real>& priorityVec) const
{
	treeVec.push_back(static_cast<int>(theFT));
	int n = children.size();
	treeVec.push_back(n);
	for(int k = 0; k < n; k++)
		children[k]->getMPIDataFormat(treeVec, ptVec, matVec, priorityVec);
}

void CSGFunctionNode::getMPIDataFormat(vector<int>& treeVec, vector<Point_3>& ptVec,
                    vector<GenericMaterial>& matVec, vector<Real>& priorityVec) const
{
    treeVec.push_back(static_cast<int>(theFT));
    int n = children.size();
    treeVec.push_back(n);
    for(int k = 0; k < n; k++)
        children[k]->getMPIDataFormat(treeVec, ptVec, matVec, priorityVec);
}

bool CSGFunctionNode::checkMeshability()
{
	int n = children.size();
	for(int k = 0; k < n; k++)
	{
		if(!children[k]->checkMeshability())
			return false;
	}
	return true;
}

Real CSGFunctionNode::getMaxPointRadius() const
{
	int n = children.size();
	Real maxRadius = 0;
	for(int k = 0; k < n; k++)
	{
		Real tmpRad = children[k]->getMaxPointRadius();
		maxRadius = MAX(tmpRad, maxRadius);
	}
	return maxRadius;
}
}

