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

#include <algorithm>
#include <iostream>
#include <unordered_set>
#include "mesh2d.h"
#include "element.h"
#include "elemedge.h"
#include "lintri.h"
#include "linquad.h"
#include "tomesh.h"

using std::vector;

Mesh2D::Mesh2D(const Topologies::TOMesh* const inMesh, const Topologies::GenericMaterial& baseMat):
	numElements(inMesh->getNumElements()),
  numEdges(0),
  numNodes(inMesh->getNumNodes()),
  numInternalEdges(0),
  numBoundaryEdges(0),
  numUnk(0)
{
	nodeVec.reserve(numNodes);
	elemVec.reserve(numElements);
	for(std::size_t k = 0; k < numNodes; ++k)
	{
		nodeVec.push_back(std::unique_ptr<Point2D>(new Point2D(inMesh->getNode2D(k))));
		basisNodeMap[nodeVec.back().get()] = k;
	}
	std::vector<double> matparams(baseMat.getNumParameters());
	for(std::size_t k = 0; k < matparams.size(); ++k)
		matparams[k] = baseMat.getParameter(k);
	for(std::size_t k = 0; k < numElements; ++k)
	{
		const std::vector<std::size_t>& tmpEV = inMesh->getElementConnectivity(k);
		std::vector<Point2D*> curPtVec(tmpEV.size());
		for(std::size_t kp = 0; kp < tmpEV.size(); ++kp)
			curPtVec[kp] = nodeVec[tmpEV[kp]].get();
		std::vector<double> curparams = matparams;
		curparams[1] *= inMesh->getOptVal(k);
		curparams[2] *= inMesh->getOptVal(k);
		Topologies::GenericMaterial curMat(curparams);
		if(tmpEV.size() == 3)
			elemVec.push_back(std::unique_ptr<Element<Point2D>>(new LinearTriangle<Point2D>(curPtVec, curMat)));
		else
			elemVec.push_back(std::unique_ptr<Element<Point2D>>(new LinearQuadrilateral<Point2D>(curPtVec, curMat)));
	}
	double area;
	bool err = cleanSurfaceDefinition(area);
	assert(!err);
}

Mesh2D::Mesh2D(bool& err, const vector<Point2D*>& inNodeVec, const vector<Element<Point2D>*>& inElemVec):
numElements(0),
	numEdges(0),
	numNodes(0),
	numInternalEdges(0),
	numBoundaryEdges(0),
	numUnk(0)
{
	finishSetup(err, inNodeVec, inElemVec);
}

void Mesh2D::finishSetup(bool& err, const vector<Point2D*>& inNodeVec, const vector<Element<Point2D>*>& inElemVec)
{
	numNodes = inNodeVec.size();
	numElements = inElemVec.size();
	nodeVec.reserve(numNodes);
	elemVec.reserve(numElements);
//	std::cout << "Copying points" << std::endl;
	for(std::size_t k = 0; k < numNodes; k++)
	{
		nodeVec.push_back(std::unique_ptr<Point2D>(new Point2D(*(inNodeVec[k]))));
		basisNodeMap[nodeVec[nodeVec.size() - 1].get()] = k;
	}
//	std::cout << "Copying elements" << std::endl;
	for(std::size_t k = 0; k < numElements; k++)
	{
		unsigned numElemNodes = inElemVec[k]->getNumNodes();
		std::vector<std::size_t> tnID(numElemNodes);
		for(unsigned n = 0; n < numElemNodes; n++)
		{
			std::vector<Point2D*>::const_iterator fit = std::find(inNodeVec.begin(), inNodeVec.end(), inElemVec[k]->getNode(n));
			assert(fit != inNodeVec.end()); // Fails if node wasn't found in inNodeVec
			tnID[n] = fit - inNodeVec.begin();
		}
		elemVec.push_back(copyElement(inElemVec[k], tnID));
	}
	double area;
	err = cleanSurfaceDefinition(area);
}

std::unique_ptr<Element<Point2D>> Mesh2D::copyElement(const Element<Point2D>* const inElem, const std::vector<std::size_t>& ptIDVec) const
{
	Topologies::GenericMaterial mat = inElem->getGenericMaterial();
	PatchType elemPT = inElem->getPatchType();
	if(elemPT == ptLinTri)
	{
		assert(ptIDVec.size() == 3);
		std::vector<Point2D*> elemConnVec(3);
		elemConnVec[0] = nodeVec[ptIDVec[0]].get();
		elemConnVec[1] = nodeVec[ptIDVec[1]].get();
		elemConnVec[2] = nodeVec[ptIDVec[2]].get();
		return std::unique_ptr<Element<Point2D>>(new LinearTriangle<Point2D>(elemConnVec, mat));
	}
	else if(elemPT == ptLinQuad)
	{
		assert(ptIDVec.size() == 4);
		std::vector<Point2D*> elemConnVec(4);
		elemConnVec[0] = nodeVec[ptIDVec[0]].get();
		elemConnVec[1] = nodeVec[ptIDVec[1]].get();
		elemConnVec[2] = nodeVec[ptIDVec[2]].get();
		elemConnVec[3] = nodeVec[ptIDVec[3]].get();
		return std::unique_ptr<Element<Point2D>>(new LinearQuadrilateral<Point2D>(elemConnVec, mat));
	}
	else
	{
		std::cout << "Error: Unrecognized element type, aborting" << std::endl;
		assert(false);
	}
	return nullptr;
}

Mesh2D::Mesh2D(const Mesh2D& copy)
{
	std::vector<Point2D*> copyNodeVec(copy.nodeVec.size());
	for(std::size_t k = 0; k < copyNodeVec.size(); k++)
		copyNodeVec[k] = copy.nodeVec[k].get();
	std::vector<Element<Point2D>*> copyElemVec(copy.elemVec.size());
	for(std::size_t k = 0; k < copyElemVec.size(); k++)
		copyElemVec[k] = copy.elemVec[k].get();
	bool err;
	finishSetup(err, copyNodeVec, copyElemVec);
}

Mesh2D::Mesh2D(Mesh2D&& copy)
{
	swap(copy);
}

Mesh2D& Mesh2D::operator=(Mesh2D rhs)
{
	swap(rhs);
	return *this;
}

void Mesh2D::swap(Mesh2D& arg2)
{
	std::swap(numElements, arg2.numElements);
	std::swap(numEdges, arg2.numEdges);
	std::swap(numNodes, arg2.numNodes);
	std::swap(numInternalEdges, arg2.numInternalEdges); 
	std::swap(numBoundaryEdges, arg2.numBoundaryEdges);
	std::swap(numUnk, arg2.numUnk);
	nodeVec.swap(arg2.nodeVec);
	edgeVec.swap(arg2.edgeVec);
	elemVec.swap(arg2.elemVec);
	basisNodeMap.swap(arg2.basisNodeMap);
}

Mesh2D::~Mesh2D()
{
}

bool Mesh2D::cleanSurfaceDefinition(double& area)
{
//	std::cout << "Setting connectivity" << std::endl;
	bool err = setConnectivity(); // set up pointers from patches to edges,
	// edges to patches, and edges to nodes.
//	std::cout << "Mesh info: num elements: " << numElements << ", num edges: " << numEdges 
//					<< ", num nodes: " << numNodes << std::endl;
	if(err)
		return err;
//	std::cout << "Sorting boundary edges" << std::endl;
	putBoundaryEdgesLast(); // Ensure that all edges on the boundary
	// are listed last in itsEdges.
//	std::cout << "Fixing Jacobians" << std::endl;
	ensureConsistentNormals(); // Ensure that all patch normals point
	// in the same direction.
	makeUpwardNormals();	// make normals point in the +z direction
//	std::cout << "Computing area of mesh ... ";
	area = computeArea();
//	std::cout << "area = " << area << std::endl;
//	std::cout << "Setting standard edge order" << std::endl;
	setStandardEdgeOrder(); // Reorder edge pointers in a standard way
//	std::cout << "Setting basis function map" << std::endl;
	setPatch2BasisMap();
	// reorderRCM();
//	outputInfo(area);
//	std::cout << "Done mesh processing" << std::endl;
	// relative to the node pointers
	return false;
}

void Mesh2D::outputInfo(double area) const
{
	std::cout << "Mesh info: ";
	std::cout << "\nNumber of Nodes: " << numNodes;
	std::cout << "\nNumber of Patches: " << numElements;
	std::cout << "\nNumber of Edges: " << numEdges;
	std::cout << "\n\tBoundary: " << numBoundaryEdges;
	std::cout << "\n\tInternal: " << numInternalEdges;
	std::cout << "\nArea (m^2): " << area;
	std::cout << "\nArea (square meters):" << area;
	std::cout << "\nNumber of total basis functions: " << numUnk;
	double bfpsw = numUnk/area;
	std::cout << "\nBasis functions per square wavelength: " << bfpsw;
	if (bfpsw < 200.)
	{
		std::cout << "\n**********************************\n";
		std::cout << "\n             WARNING!             \n";
		std::cout << "\n**********************************\n";
		std::cout << "\nFewer than 200 bases per square wavelength.";
	}
	std::cout << std::endl << std::endl;
}	

double Mesh2D::computeArea() const
{
	double area = 0.;
	for (std::size_t iTri = 0; iTri < numElements; iTri++)
		area += elemVec[iTri]->getArea();
	return area;
}

double Mesh2D::computeMaterialPropertyTotal(unsigned propID) const
{
	double propTot = 0.;
	for (std::size_t iTri = 0; iTri < numElements; iTri++)
	{
		Topologies::GenericMaterial tmp = elemVec[iTri]->getGenericMaterial();
		propTot += elemVec[iTri]->getArea()*tmp.getParameter(propID);
	}

	return propTot;
}

bool Mesh2D::setConnectivity()
{
	edgeVec.clear();
	edgeVec.reserve(10*numElements);
	std::unique_ptr<ElemEdge<Point2D>> nextEdge;
	bool err = true;
	typedef std::unordered_set<ElemEdge<Point2D>*, Elemedge_hash<Point2D>, Elemedge_equal<Point2D>> Elemedge_set;
	Elemedge_set eset;
	for(std::size_t iElem = 0; iElem < numElements; iElem++)
	{
		for(unsigned short int jEdge = 0; jEdge < elemVec[iElem]->getNumNodes(); jEdge++)
		{
			nextEdge = elemVec[iElem]->createEdge(jEdge);
			Elemedge_set::const_iterator eit = eset.find(nextEdge.get()); // Search for edge
			if(eit == eset.end())
			{
				err = nextEdge->addElement(elemVec[iElem].get());
				eset.insert(nextEdge.get());
				edgeVec.push_back(std::move(nextEdge));
			}
			else
				err = (*eit)->addElement(elemVec[iElem].get());
			if(err)
				return err;
		} // End loop over local edges of the patch
	} // End loop over patches
	edgeVec.shrink_to_fit();
	numEdges = edgeVec.size();
	// Add edge list to data structure, and create patch to edge pointers
	Element<Point2D>* pElem = nullptr;
	ElemEdge<Point2D>* theEdge = nullptr;
	for (std::size_t i = 0; i < numEdges; i++)
	{
		theEdge = edgeVec[i].get();
		pElem = edgeVec[i]->elem0();
		pElem->addEdge(theEdge);
		pElem = edgeVec[i]->elem1();
		if(pElem != nullptr) // Boundary edges don't have a second element
			pElem->addEdge(theEdge);
	}
	return err;
}

void Mesh2D::putBoundaryEdgesLast()
{
	numBoundaryEdges = 0;
	numInternalEdges = 0;
	do
	{
		if(edgeVec[numInternalEdges]->isBoundaryEdge())
		{
			++numBoundaryEdges; // Increment number of boundary edges
			edgeVec[numInternalEdges].swap(edgeVec[numEdges - numBoundaryEdges]);
		}
		else
			++numInternalEdges;
	} while(numBoundaryEdges + numInternalEdges < numEdges);
}

Point2D* Mesh2D::findBoundaryPoint(Point2D pt)
{
	for(std::size_t k = numInternalEdges + 1; k < numEdges; k++)
	{
		ElemEdge<Point2D>* tmpTE = edgeVec[k].get();
		Point2D tmpPt0 = *(tmpTE->node0());
		if(pt == tmpPt0)
			return tmpTE->node0();
		tmpPt0 = *(tmpTE->node1());
		if(pt == tmpPt0)
			return tmpTE->node1();
	}
	return 0;
}

void Mesh2D::ensureConsistentNormals()
{
	for(std::size_t kelem = 0; kelem < numElements; kelem++)
	{
		std::stack<Element<Point2D>*> theStack;
		if(!elemVec[kelem]->isReordered())
		{
			elemVec[kelem]->hasBeenReordered();
			theStack.push(elemVec[kelem].get());

			Element<Point2D>* pParent = nullptr;
			Element<Point2D>* pChild = nullptr;
			ElemEdge<Point2D>* pParentEdge = nullptr;
			std::unique_ptr<ElemEdge<Point2D>> ParentEdge, ChildEdge;
			unsigned short int jEdge, kEdge, temp;
			temp = 0;

			while(!theStack.empty())
			{
				pParent = theStack.top();
				theStack.pop();
				for(jEdge = 0; jEdge < pParent->getNumNodes(); jEdge++)
				{
					pParentEdge = pParent->getEdge(jEdge);
					for(kEdge = 0; kEdge < pParent->getNumNodes(); kEdge++)
					{
						ParentEdge = pParent->createEdge(kEdge);
						if(*ParentEdge == *pParentEdge)
							break;
					}	
					pChild = pParentEdge->elem0();
					if(pChild == pParent)
						pChild = pParentEdge->elem1();
					if(pChild == 0)
						continue;
					if(pChild->isReordered())
						continue;
					for(kEdge = 0; kEdge < pChild->getNumNodes(); kEdge++)
					{
						ChildEdge = pChild->createEdge(kEdge);
						if(*ChildEdge == *ParentEdge)
							break;
					}
					if(!ParentEdge->isAlignedWith(*ChildEdge))
						pChild->switchXiEta();

					pChild->hasBeenReordered();
					theStack.push(pChild);
				}
			}
		}
	}     
}  

void Mesh2D::makeUpwardNormals()
{
	for(std::size_t iElem = 0; iElem < numElements; iElem++)
	{
		if(elemVec[iElem]->getOrientation() < 0)
			elemVec[iElem]->switchXiEta();
	}
}

void Mesh2D::setStandardEdgeOrder()
{
	for(std::size_t iElem = 0; iElem < numElements; iElem++)
		elemVec[iElem]->setStandardEdgeOrder();
}

void Mesh2D::setPatch2BasisMap()
{
	numUnk = numNodes;
	setNodeBFMap();
}

void Mesh2D::setNodeBFMap()
{
	for(std::size_t iElem = 0; iElem < numElements; iElem++)
	{
		Element<Point2D>* curElem = elemVec[iElem].get();
		for(std::size_t jNode = 0; jNode < curElem->getNumNodes(); ++jNode)
		{
			Point2D* curNode = curElem->getNode(jNode);
			std::size_t bfid = basisNodeMap[curNode];
			curElem->addNodeBases(curNode, bfid);
		}
	}
}

std::size_t Mesh2D::findNode(const Point2D* theNode) const
{
	bool found = false;
	std::size_t nodeID = 0;
	for(std::size_t k = 0; k < numNodes && !found; k++)
	{
		if(nodeVec[k].get() == theNode)
		{
			nodeID = k;
			found = true;
		}
	}
	assert(found);
	return nodeID;
}

void Mesh2D::reorderRCM()
{
	// Apply reverse Cuthill-McKee reordering to ensure a banded matrix
	bool *nodeMarkList = new bool[numNodes];
	int *permIDList = new int[numNodes];
	std::size_t kount = 0;
	for(std::size_t k = 0; k < numNodes; k++)
		nodeMarkList[k] = false;
	// put all boundary nodes first:
	std::list<Point2D*> curLevelSet;
	std::list<int> curLevelSetDegrees;
	Point2D* tmpPt1 = edgeVec[numInternalEdges]->node0();
	int tmpID1 = findNode(tmpPt1);
	insertSortDegree(tmpPt1, curLevelSet, curLevelSetDegrees);
	permIDList[kount] = tmpID1;
	nodeMarkList[tmpID1] = true;
	kount++;

	tmpPt1 = edgeVec[numInternalEdges]->node1();
	tmpID1 = findNode(tmpPt1);
	insertSortDegree(tmpPt1, curLevelSet, curLevelSetDegrees);
	permIDList[kount] = tmpID1;
	nodeMarkList[tmpID1] = true;
	kount++;

	for(std::size_t k = 1; k < numBoundaryEdges; k++)
	{
		Point2D *pt0 = edgeVec[numInternalEdges + k]->node0(), *pt1 = edgeVec[numInternalEdges + k]->node1();
		int ptID0 = findNode(pt0), ptID1 = findNode(pt1);
		if(!nodeMarkList[ptID0])
		{
			tmpPt1 = pt0;
			tmpID1 = findNode(tmpPt1);
			insertSortDegree(tmpPt1, curLevelSet, curLevelSetDegrees);
			permIDList[kount] = tmpID1;
			nodeMarkList[tmpID1] = true;
			kount++;
		}
		if(!nodeMarkList[ptID1])
		{
			tmpPt1 = pt1;
			tmpID1 = findNode(tmpPt1);
			insertSortDegree(tmpPt1, curLevelSet, curLevelSetDegrees);
			permIDList[kount] = tmpID1;
			nodeMarkList[tmpID1] = true;
			kount++;
		}
	}

	// traverse graph breadth-first
	std::size_t nodeCount = 0;
	while(kount < numNodes && nodeCount < 2*numNodes)
	{
		std::list<Point2D*> nextLevelSet;
		std::list<int> nextLevelSetDegrees;
		std::list<Point2D*>::iterator pPtList = curLevelSet.begin();
		for(; pPtList != curLevelSet.end(); ++pPtList)
		{
			vector<ElemEdge<Point2D>*> edgeVec;
			Point2D* curPt = *pPtList;
			getNeighboringEdges(curPt, edgeVec);
			for(std::size_t k = 0; k < edgeVec.size(); k++)
			{
				Point2D* neighborPt;
				if(edgeVec[k]->node0() == curPt)
					neighborPt = edgeVec[k]->node1();
				else
					neighborPt = edgeVec[k]->node0();

				int neighborID = findNode(neighborPt);
				if(!nodeMarkList[neighborID])
				{
					nodeMarkList[neighborID] = true;
					insertSortDegree(neighborPt, nextLevelSet, nextLevelSetDegrees);
					permIDList[kount] = neighborID;
					kount++;
				}
			}
		}
		curLevelSet = nextLevelSet;
		curLevelSetDegrees = nextLevelSetDegrees;
		nodeCount++;
	}
	if(kount != numNodes)
	{
		std::cout << "error in reordering!" << std::endl;
		for(std::size_t k = 0; k < numNodes; k++)
		{
			if(!nodeMarkList[k])
			{
				permIDList[k] = k;
				nodeMarkList[k] = true;
			}
		}
	}
	// reorder nodes:
	std::vector<std::unique_ptr<Point2D>> oldNodeList;
	oldNodeList.reserve(nodeVec.size());
	for(std::size_t k = 0; k < numNodes; k++)
		oldNodeList[k] = std::move(nodeVec[k]);
	for(std::size_t k = 0; k < numNodes; k++)
		nodeVec[k] = std::move(oldNodeList[permIDList[numNodes - k - 1]]);
	// remap bf ids:
	setNodeBFMap();
	// Clean up
	delete [] nodeMarkList;
	delete [] permIDList;
}

int Mesh2D::getNodeDegree(Point2D* tstPt)
{
	int degree = 0;
	for(std::size_t k = 0; k < numEdges; k++)
	{
		if(edgeVec[k]->node0() == tstPt || edgeVec[k]->node1() == tstPt)
			degree++;
	}
	return degree;
}

void Mesh2D::getNeighboringEdges(Point2D* tstPt, vector<ElemEdge<Point2D>*>& edgeVec)
{
	for(std::size_t k = 0; k < numEdges; k++)
	{
		if(edgeVec[k]->node0() == tstPt || edgeVec[k]->node1() == tstPt)
			edgeVec.push_back(edgeVec[k]);
	}
}

void Mesh2D::insertSortDegree(Point2D* inPt, std::list<Point2D*>& levelSet, 
	std::list<int>& levelSetDegrees)
{
	int newDegree = getNodeDegree(inPt);
	std::list<Point2D*>::iterator pPtList = levelSet.begin();
	std::list<int>::iterator pDegreeList = levelSetDegrees.begin();
	while(pDegreeList != levelSetDegrees.end() && *pDegreeList > newDegree)
	{
		++pDegreeList;
		++pPtList;
	}
	levelSet.insert(pPtList, inPt);
	levelSetDegrees.insert(pDegreeList, newDegree);
}

void Mesh2D::printMesh(const std::string& fileName) const
{
	std::ofstream outFile(fileName.c_str());
	outFile << "nodes = [";
	for(std::size_t k = 0; k < numNodes; k++)
		outFile << nodeVec[k]->x << " " << nodeVec[k]->y << std::endl;
	outFile << "];" << std::endl;
	outFile << "elems = [";
	for(std::size_t k = 0; k < numElements; k++)
	{
		Element<Point2D>* curElem = elemVec[k].get();
		for(unsigned kc = 0; kc < curElem->getNumNodes(); kc++)
		{
			Point2D* cNode = curElem->getNode(kc);
			outFile << findNode(cNode)+1 << " ";
		}
		outFile << std::endl;
	}
	outFile << "];" << std::endl;
	outFile << "mats = [";
	for(std::size_t k = 0; k < numElements; k++)
	{
		Topologies::GenericMaterial mat = elemVec[k]->getGenericMaterial();
		int nm = mat.getNumParameters();
		for(int k2 = 0; k2 < nm; k2++)
			outFile << mat.getParameter(k2) << " ";

		outFile << std::endl;
	}
	outFile << "];" << std::endl;
}

void Mesh2D::printBoundaries(std::ofstream& outFile) const
{
	outFile << "hold on;" << std::endl;
	for(std::size_t k = 0; k < numEdges; k++)
	{
		Element<Point2D>* tmpTri0 = edgeVec[k]->elem0();
		Element<Point2D>* tmpTri1 = edgeVec[k]->elem1();
		Point2D p0 = *(edgeVec[k]->node0()), p1 = *(edgeVec[k]->node1());
		if(tmpTri1 == nullptr)
		{
			outFile << "plot([" << p0.x << " " << p1.x << "], [" << p0.y << " " << p1.y 
				<< "], 'k', 'LineWidth', 2);" << std::endl;
		}
		else
		{
			double epsR0 = tmpTri0->getGenericMaterial().getParameter(0), epsR1 = tmpTri1->getGenericMaterial().getParameter(0);
			if(epsR0 != epsR1)
			{
				outFile << "plot([" << p0.x << " " << p1.x << "], [" << p0.y << " " << p1.y 
					<< "], 'k', 'LineWidth', 2);" << std::endl;
			}
		}
	}
}
