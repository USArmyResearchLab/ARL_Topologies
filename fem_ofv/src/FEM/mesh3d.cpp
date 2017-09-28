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
#include <unordered_map>
#include <unordered_set>
#include "mesh3d.h"
#include "element.h"
#include "point3d.h"
#include "cell.h"
#include "lintetra.h"
#include "trilinhex.h"
#include "tomesh.h"

Mesh3D::Mesh3D(const Topologies::TOMesh* const inMesh, const Topologies::GenericMaterial& baseMat)
{
	itsNumNodes = inMesh->getNumNodes();
	itsNumCells = inMesh->getNumElements();
	nodeVec.reserve(itsNumNodes);
	cellVec.reserve(itsNumCells);
	// Generate list of all grid points
	for(std::size_t k = 0; k < itsNumNodes; ++k)
	{
		nodeVec.push_back(std::unique_ptr<Point3D>(new Point3D(inMesh->getNode3D(k))));
		basisNodeMap[nodeVec.back().get()] = k;
	}
	// Set up FEM mesh
	std::vector<double> matparams(baseMat.getNumParameters());
	for(std::size_t k = 0; k < matparams.size(); ++k)
		matparams[k] = baseMat.getParameter(k);
	for(std::size_t k = 0; k < itsNumCells; ++k)
	{
		const std::vector<std::size_t>& tmpEV = inMesh->getElementConnectivity(k);
		std::vector<Point3D*> curElemVec(tmpEV.size());
		for(std::size_t kp = 0; kp < tmpEV.size(); ++kp)
			curElemVec[kp] = nodeVec[tmpEV[kp]].get();
		std::vector<double> curparams = matparams;
		curparams[1] *= inMesh->getOptVal(k);
		curparams[2] *= inMesh->getOptVal(k);
		Topologies::GenericMaterial curMat(curparams);
		if(tmpEV.size() == 4)
			cellVec.push_back(std::unique_ptr<Cell>(new LinTetra(curElemVec, curMat)));
		else
			cellVec.push_back(std::unique_ptr<Cell>(new TriLinHex(curElemVec, curMat)));
	}
	double area, volume;
	// Sets up the cell connectivity, patch list, basis function mapping
	cleanVolumeDefinition(area, volume);
}

Mesh3D::Mesh3D(bool& err, const std::vector<Point3D*>& inNodeVec, const std::vector<Cell*>& inElemVec)
{
	finishSetup(err, inNodeVec, inElemVec);
}

void Mesh3D::finishSetup(bool& err, const std::vector<Point3D*>& inNodeVec, const std::vector<Cell*>& inElemVec)
{
	itsNumNodes = inNodeVec.size();
	itsNumCells = inElemVec.size();
	nodeVec.reserve(itsNumNodes);
	cellVec.reserve(itsNumCells);
	// Copy nodes
//	std::cout << "Copying nodes." << std::endl;
	std::unordered_map<Point3D*, std::size_t> pointPtrMap; // Used for fast searching
	for(std::size_t k = 0; k < itsNumNodes; k++)
	{
		nodeVec.push_back(std::unique_ptr<Point3D>(new Point3D(*(inNodeVec[k]))));
		pointPtrMap[inNodeVec[k]] = k;
		basisNodeMap[nodeVec[nodeVec.size() - 1].get()] = k;
	}
	// Copy elements
//	std::cout << "Copying " << itsNumCells << " elements." << std::endl;
	for(std::size_t k = 0; k < itsNumCells; k++)
	{
		std::size_t numElemNodes = inElemVec[k]->getNumNodes();
		std::vector<std::size_t> tnID(numElemNodes);
		for(std::size_t n = 0; n < numElemNodes; n++)
		{
			std::unordered_map<Point3D*, std::size_t>::const_iterator fit = pointPtrMap.find(inElemVec[k]->getNode(n));
			assert(fit != pointPtrMap.end()); // Fails if node wasn't found in inNodeVec
			tnID[n] = fit->second;
		}
		cellVec.push_back(copyCell(inElemVec[k], tnID));
	}
	double area, volume;
	// Sets up the cell connectivity, patch list, basis function mapping
	cleanVolumeDefinition(area, volume);
//	outputInfo(area, volume);
}

Mesh3D::Mesh3D(const Mesh3D& copy)
{
	std::vector<Point3D*> copyNodeVec(copy.nodeVec.size());
	for(std::size_t k = 0; k < copyNodeVec.size(); k++)
		copyNodeVec[k] = copy.nodeVec[k].get();
	std::vector<Cell*> copyElemVec(copy.cellVec.size());
	for(std::size_t k = 0; k < copyElemVec.size(); k++)
		copyElemVec[k] = copy.cellVec[k].get();
	bool err;
	finishSetup(err, copyNodeVec, copyElemVec);
}

Mesh3D::Mesh3D(Mesh3D&& copy)
{
	swap(copy);
}

Mesh3D& Mesh3D::operator=(Mesh3D rhs)
{
	swap(rhs);
	return *this;
}

void Mesh3D::swap(Mesh3D& arg2)
{
	std::swap(itsNumCells, arg2.itsNumCells);
	std::swap(itsNumPatches, arg2.itsNumPatches);
	std::swap(itsNumEdges, arg2.itsNumEdges);
	std::swap(itsNumNodes, arg2.itsNumNodes);
	std::swap(itsFaceNunk, arg2.itsFaceNunk);
	std::swap(itsInternalNunk, arg2.itsInternalNunk);
	std::swap(typeOfCell, arg2.typeOfCell);
	std::swap(itsNumInternalPatches, arg2.itsNumInternalPatches);
	std::swap(itsNumBoundaryPatches, arg2.itsNumBoundaryPatches);
	std::swap(itsNunk, arg2.itsNunk);
	std::swap(itsOffset, arg2.itsOffset);
	nodeVec.swap(arg2.nodeVec);
	edgeVec.swap(arg2.edgeVec);
	patchVec.swap(arg2.patchVec);
	cellVec.swap(arg2.cellVec);
	basisNodeMap.swap(arg2.basisNodeMap);
}

Mesh3D::~Mesh3D()
{
}

void Mesh3D::setOffset(std::size_t off)
{
	assert(itsOffset == 0);
	itsOffset = off;
}

std::unique_ptr<Cell> Mesh3D::copyCell(const Cell* const inCell, const std::vector<std::size_t>& ptIDVec) const
{
	Topologies::GenericMaterial mat = inCell->getGenericMaterial();
	CellType cellCT = inCell->getCellType();
	if(cellCT == ctLinTetra)
	{
		assert(ptIDVec.size() == 4);
		std::vector<Point3D*> elemConnVec(4);
		elemConnVec[0] = nodeVec[ptIDVec[0]].get();
		elemConnVec[1] = nodeVec[ptIDVec[1]].get();
		elemConnVec[2] = nodeVec[ptIDVec[2]].get();
		elemConnVec[3] = nodeVec[ptIDVec[3]].get();
		return std::unique_ptr<Cell>(new LinTetra(elemConnVec, mat));
	}
	else if(cellCT == ctTriLinHex)
	{
		assert(ptIDVec.size() == 8);
		std::vector<Point3D*> elemConnVec(8);
		elemConnVec[0] = nodeVec[ptIDVec[0]].get();
		elemConnVec[1] = nodeVec[ptIDVec[1]].get();
		elemConnVec[2] = nodeVec[ptIDVec[2]].get();
		elemConnVec[3] = nodeVec[ptIDVec[3]].get();
		elemConnVec[4] = nodeVec[ptIDVec[4]].get();
		elemConnVec[5] = nodeVec[ptIDVec[5]].get();
		elemConnVec[6] = nodeVec[ptIDVec[6]].get();
		elemConnVec[7] = nodeVec[ptIDVec[7]].get();
		return std::unique_ptr<Cell>(new TriLinHex(elemConnVec, mat));
	}
	std::cout << "Error: Unrecognized element type, aborting" << std::endl;
	assert(false);
	return nullptr;
}

void Mesh3D::cleanVolumeDefinition(double& area, double& volume)
{
//	std::cout << "Setting patch connectivity" << std::endl;
	setPatchConnectivity(); // set up pointers from cell to patches, patches to cells, 
	//patches to nodes.
//	std::cout << "Setting edge connectivity" << std::endl;
	setEdgeConnectivity(); // set up pointers from patches to edges, etc.
//	std::cout << "Nodes: " << itsNumNodes << "\nEdges: " << itsNumEdges << "\nFaces: " << itsNumPatches << "\nCells: " << itsNumCells << std::endl;
//	std::cout << "Euler number: " << itsNumNodes - itsNumEdges + itsNumPatches - itsNumCells << std::endl;
//	std::cout << "Sorting boundary patches" << std::endl;
	putBoundaryPatchesLast(); // Ensure that all patches on the boundary
	// are listed last in itspatches.
//	std::cout << "Setting jacobian signs" << std::endl;
	ensurePositiveJacobian(); // Ensure that the Jacobian of each cell is positive
	area = computeArea(); // computes the surface area by just looping over the bounadry patches
	volume = computeVolume(); // Volume computation 
//	std::cout << "Area: " << area << ", volume: " << volume << std::endl;
//	std::cout << "Fixing node orderings" << std::endl;
	setStandardPatchOrder(); // Reorder patch pointers in a standard way
	// relative to the node pointers
//	std::cout << "Setting basis function numbers" << std::endl;
	setCell2BasisMap(); // maps the face and internal basis functions to the each cell 
}

void Mesh3D::ensurePositiveJacobian()
{
	double tempVol;
	for (std::size_t iCell = 0; iCell < itsNumCells; iCell++)
	{
		tempVol = cellVec[iCell]->volumeIntegral();
		if (tempVol < 0.)
			cellVec[iCell]->switchXi2Xi3();
	}	
}

double Mesh3D::computeArea() const
{
	double area = 0.;
	for (std::size_t iPatch = itsNumInternalPatches; iPatch < itsNumPatches; iPatch++)
		area += fabs(patchVec[iPatch]->getArea());
	return area;
}

double Mesh3D::computeVolume() const
{
	double volume = 0.;
	for (std::size_t iCell = 0; iCell < itsNumCells; iCell++)
		volume += cellVec[iCell]->volumeIntegral();
	return volume;
}

bool Mesh3D::setPatchConnectivity()
{
	// Creates a list of patches based on the cells that make up the volume. Then the   
	// patches are added to the cell structures.  In the end, each patch structure 
	// knows which 2 cells it joins and each cell knows each of its patches.
	itsNumPatches = 0;
	patchVec.clear();
	patchVec.reserve(6*itsNumCells);
	std::unique_ptr<Element<Point3D>> nextPatch;
	typedef std::unordered_set<Element<Point3D>*, Element_hash<Point3D>, Element_equal<Point3D>> Element_set;
	Element_set pset; // Hash table for fast searching of existing patches
	for (std::size_t iCell = 0; iCell < itsNumCells; iCell++)
	{
		for (unsigned short int jPatch = 0; jPatch < cellVec[iCell]->getNumPatches(); jPatch++)
		{	
			// creates a patch using the 3 cell nodes 
			nextPatch = cellVec[iCell]->createPatch(jPatch);	
			Element_set::const_iterator pit = pset.find(nextPatch.get()); // Search for patch
			if(pit == pset.end())
			{
				if(nextPatch->addCell(cellVec[iCell].get()))
				{
					std::cout << "Saw error in Element::addCell" << std::endl;
					return true;
				}
				pset.insert(nextPatch.get()); // Add patch into hash table for searching
				patchVec.push_back(std::move(nextPatch)); // Store unique_ptr in patchVec
			}
			else
				(*pit)->addCell(cellVec[iCell].get());
		} // End loop over local patches of the cell
	} // End loop over cells
	patchVec.shrink_to_fit();
	itsNumPatches = patchVec.size();
	// Add patch list to data structure, and create cell to patch pointers
	Cell* pCell = nullptr;
	Element<Point3D>* pPatch = nullptr;
	for (std::size_t i = 0; i < itsNumPatches; i++)
	{
		pPatch = patchVec[i].get();
		pCell = pPatch->cell0();		
		// adding the patch details to each of the cell it belongs to
		pCell->addPatch(pPatch);
		pCell = pPatch->cell1();
		if (pCell != 0)
			pCell->addPatch(pPatch);
	}
	return false;
}

bool Mesh3D::setEdgeConnectivity()
{
	edgeVec.clear();
	edgeVec.reserve(4*itsNumPatches);
	std::unique_ptr<ElemEdge<Point3D>> nextEdge;
	bool err = false;
	typedef std::unordered_set<ElemEdge<Point3D>*, Elemedge_hash<Point3D>, Elemedge_equal<Point3D>> Elemedge_set;
	Elemedge_set eset;
	for(std::size_t iElem = 0; iElem < itsNumPatches; iElem++)
	{
		for(unsigned short int jEdge = 0; jEdge < patchVec[iElem]->getNumNodes(); jEdge++)
		{
			nextEdge = patchVec[iElem]->createEdge(jEdge);
			Elemedge_set::const_iterator eit = eset.find(nextEdge.get()); // Search for edge
			if(eit == eset.end())
			{
				nextEdge->addElement(patchVec[iElem].get());
				eset.insert(nextEdge.get());
				edgeVec.push_back(std::move(nextEdge));
			}
			else
				(*eit)->addElement(patchVec[iElem].get());
		} // End loop over local edges of the patch
	} // End loop over patches
	edgeVec.shrink_to_fit();
	itsNumEdges = edgeVec.size();
	// Add edge list to data structure, and create patch to edge pointers
	Element<Point3D>* pElem = nullptr;
	ElemEdge<Point3D>* theEdge = nullptr;
	for (std::size_t i = 0; i < itsNumEdges; i++)
	{
		theEdge = edgeVec[i].get();
		for(unsigned j = 0; j < theEdge->getNumAttachedElems(); ++j)
		{
			pElem = edgeVec[i]->elemN(j);
			if(pElem != nullptr)
				pElem->addEdge(theEdge);
		}
	}
	return err;
}

void Mesh3D::putBoundaryPatchesLast()
{
	//Sorts the array of patches so that the boundary patches are listed last.
	itsNumBoundaryPatches = 0;
	itsNumInternalPatches = 0;
	do
	{
		if (patchVec[itsNumInternalPatches]->isBoundaryPatch())
		{
			++itsNumBoundaryPatches; // Increment number of boundary patches
			patchVec[itsNumInternalPatches].swap(patchVec[itsNumPatches - itsNumBoundaryPatches]);
		}
		else
			++itsNumInternalPatches;
	} while (itsNumBoundaryPatches + itsNumInternalPatches < itsNumPatches);
}

void Mesh3D::setStandardPatchOrder()
{
	for(std::size_t iCell = 0; iCell < itsNumCells; iCell++)
		cellVec[iCell]->setStandardPatchOrder();
}

void Mesh3D::setCell2BasisMap()
{
	// set up the basis functions
	itsNunk = itsNumNodes;
	setNodeBFMap();
	// Later add patch bfs and cell bfs if necessary for higher order elements
	itsFaceNunk = 0;
	itsInternalNunk = 0;
}

void Mesh3D::setNodeBFMap()
{
	for(std::size_t iElem = 0; iElem < itsNumCells; iElem++)
	{
		Cell* curElem = cellVec[iElem].get();
		for(unsigned jNode = 0; jNode < curElem->getNumNodes(); ++jNode)
		{
			Point3D* curNode = curElem->getNode(jNode);
			std::size_t bfid = basisNodeMap[curNode];
			curElem->addNodeBases(curNode, bfid);
		}
	}
}

void Mesh3D::setPatchBFMap()
{
	// Currently not implemented, only linear tetrahedra are, which don't have internal bfs
}

void Mesh3D::setCellBFMap()
{
	// Currently not implemented, only linear tetrahedra are, which don't have internal bfs
}

std::size_t Mesh3D::findCell(const Cell* theCell) const
{
	for (std::size_t i = 0; i < itsNumCells; ++i)
		if (theCell == cellVec[i].get())
			return i;
	assert(false);
	return itsNumCells;
}

std::size_t Mesh3D::findNode(const Point3D* theNode) const
{
	bool found = false;
	std::size_t nodeID = 0;
  for(std::size_t k = 0; k < itsNumNodes && !found; k++)
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

void Mesh3D::outputInfo(double area, double volume) const
{
	std::cout << "\nNumber of Nodes: " << itsNumNodes;
	std::cout << "\nNumber of Cells: " << itsNumCells;
	std::cout << "\nNumber of Patches: " << itsNumPatches;
	std::cout << "\n\tBoundary: " << itsNumBoundaryPatches;
	std::cout << "\n\tInternal: " << itsNumInternalPatches;
	std::cout << "\nNumber of total basis functions: " << itsNunk;

	std::cout << "\n\tface basis functions: " << itsFaceNunk;
	std::cout << "\n\tinternal basis functions: " << itsInternalNunk;
	std::cout << "\nvolume (cubic meters) = " << fabs(volume);
	std::cout << "\narea (square meters) = " << fabs(area);
	std::cout << "\nvolume:" << volume;
	std::cout << "\nNumber of total basis functions: " << itsNunk;
	double bfpsw = itsNunk/volume;
	std::cout << "\nBasis functions per cubic meter: " << bfpsw;
	if (bfpsw < 3000.)
	{
		std::cout << "\n**********************************\n";
		std::cout << "\n             WARNING!             \n";
		std::cout << "\n**********************************\n";
		std::cout << "\nFewer than 3000 bases per cubic wavelength.";
	}
	std::cout << std::endl << std::endl;
}	

void Mesh3D::printMesh(const std::string& fileName) const
{
	std::ofstream outFile(fileName.c_str());
	outFile << "nodes = [";
	for(std::size_t k = 0; k < itsNumNodes; k++)
		outFile << nodeVec[k]->x << " " << nodeVec[k]->y << " " << nodeVec[k]->z << std::endl;
	outFile << "];" << std::endl;
	outFile << "elems = [";
	for(std::size_t k = 0; k < itsNumCells; ++k)
	{
		Cell* curElem = cellVec[k].get();
		for(unsigned kc = 0; kc < curElem->getNumNodes(); ++kc)
		{
			Point3D* curNode = curElem->getNode(kc);
			outFile << findNode(curNode) + 1 << " " ;
		}
		outFile << std::endl;
	}
	outFile << "];" << std::endl;
	outFile << "mats = [";
	for(std::size_t k = 0; k < itsNumCells; ++k)
	{
		Topologies::GenericMaterial mat = cellVec[k]->getGenericMaterial();
		std::size_t nm = mat.getNumParameters();
		for(std::size_t k2 = 0; k2 < nm; ++k2)
			outFile << mat.getParameter(k2) << " ";
		outFile << std::endl;
	}
	outFile << "];" << std::endl;
	outFile.close();
}
