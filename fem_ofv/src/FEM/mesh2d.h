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

#ifndef MESH2D_H
#define MESH2D_H

#include <fstream>
#include <string>
#include <stack>
#include <map>
#include <unordered_map>
#include <vector>
#include <list>
#include <memory>
#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "point2d.h"
#include "UTIL/genericmaterial.h"
#include "elemedge.h"
#include "element.h"
#include "femmesh.h"

namespace Topologies{
class TOMesh;
}

//! A 2d surface mesh object that completely describes one connected surface (2D area)
/*!
Mesh2D contains all elements (which are made up of edges and nodes) and connectivity information. It can also access element matrix construction functions that are defined in the element types.  This class basically contains all necessary geometric, topological, and numerical information to describe a finite element problem.
*/
class Mesh2D : public FEMMesh
{
public:
	//! @name Constructors and destructor
	//@{
	//! Constructor which takes a TOMesh and copies the nodes and connectivity from it
	Mesh2D(const Topologies::TOMesh* const inMesh, const Topologies::GenericMaterial& baseMat);
	//! Constructor which takes in a vector of Point2D pointers and a vector of Element objects. err is returned with any error conditions
	/*! This constructor generates a 2D mesh from a set of nodes (inNodeVec) and elements (inElemVec).  These objects will be copied into for internal use and the mesh will be processed to generate and store topological information such as element edges.  err indicates whether an error has occured on return during the mesh processing.
	*/
	Mesh2D(bool& err, const std::vector<Point2D*>& inNodeVec, const std::vector<Element<Point2D>*>& inElemVec);
	Mesh2D(const Mesh2D&);
	Mesh2D(Mesh2D&&);
	Mesh2D& operator=(Mesh2D);
	void swap(Mesh2D&);
	virtual ~Mesh2D();
	//@}
	//! @name FEM element matrix functions
	//@{
	//! Returns the plane stress, linear elastic element matrix for element k
	virtual Eigen::MatrixXd getElementMatrix(std::size_t k) const {return getElement(k)->getElemMat(eptPlaneStress);}
	//! Returns either the plane strain or plane stress, linear elastic element matrix for element k
	virtual Eigen::MatrixXd getElementMatrix(std::size_t k, const ElasticProblemType theEPT) const {return getElement(k)->getElemMat(theEPT);}
	//! Returns the Laplacian element matrix for given material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(std::size_t k, double matVal) const {return getElement(k)->getLaplacianElemMat(matVal);}
	//! Returns the Laplacian element matrix for given anisotropic material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(std::size_t k, const Eigen::MatrixXd& matVal) const {return getElement(k)->getLaplacianElemMat(matVal);}
	//! Returns the stiffness matrix associated with isotropic thermal expansion
	virtual Eigen::MatrixXd getThermalExpansionMat(std::size_t k, double alphaTE, const ElasticProblemType inEPT = eptPlaneStress) const {return getElement(k)->getThermalExpansionMat(inEPT, alphaTE);}
	//! Returns a vector containing the basis function values at the given local are coordinate (xi, eta, zeta)
	virtual Eigen::VectorXd getInterpolationVec(std::size_t k, double xi, double eta, double zeta = 0.) const {return getElement(k)->getInterpolationVec(xi, eta);}
	//! Returns a matrix containing the gradient of the basis functions at element k, computed at the centroid
	virtual Eigen::MatrixXd getGradBF(std::size_t k) const {return getElement(k)->getGradBF();}
	//! Returns a matrix giving the Laplacian of the basis functions, not implemented in 3d
	virtual Eigen::MatrixXd getLaplacianBF1(std::size_t k) const {return getElement(k)->getLaplacianBF1();}
	//! Returns a matrix giving the Laplacian of the basis functions, not implemented in 3d
	virtual Eigen::MatrixXd getLaplacianBF2(std::size_t k) const {return getElement(k)->getLaplacianBF2();}
	//! Returns the mass matrix
	/*! The mass matrix is the inner product of the basis & testing functions over an element, integrated with the density.
	*/
	virtual Eigen::MatrixXd getMassMatrix(std::size_t k) const {return getElement(k)->getMassMat();}
	//! Returns the integral of the testing functions over element k
	virtual Eigen::MatrixXd getRHSTestingMatrix(std::size_t k) const {return getElement(k)->getTFInt();}
	//@}
	//! Outputs the mesh to file fileName using a Matlab-like format
	virtual void printMesh(const std::string& fileName) const;
	//! Copies the mesh 
	virtual std::unique_ptr<FEMMesh> clone() const {return std::unique_ptr<FEMMesh>(new Mesh2D(*this));}
	//! Returns the dimension of the mesn (2)
	virtual unsigned getDim() const {return 2;}
	//! Returns the material properties at element k
	virtual Topologies::GenericMaterial getGenericMaterial(std::size_t k) const {return getElement(k)->getGenericMaterial();}

	//!@name Data access, accessible from abstract class
	//@{
	//! Returns the number of elements (Cell objects) in the mesh
	virtual std::size_t getNumElements() const {return numElements;}
	//! Returns the number of nodes (or points)
	virtual std::size_t getNumUnknowns() const {return numUnk;}
	//! Returns the global basis function id associated with element kelem for local basis function klocbf
	virtual std::size_t getGlobalBF(std::size_t kelem, unsigned klocbf) const {return getElement(kelem)->getGlobalBF(klocbf);}
	//! Returns the number of nodes associated with element k
	virtual unsigned getNumElementNodes(std::size_t k) const {return getElement(k)->getNumNodes();}
	//!@name Data access, accesible only from Mesh2D
	//@{
	//! Returns the number of mesh edges
	std::size_t getNumEdges() const;
	//! Returns the number of nodes (points) in the mesh
	std::size_t getNumNodes() const;
	//! Returns the number of non-boundary edges
	std::size_t getNumInternalEdges() const;
	//! Returns the number of boundary edges
	std::size_t getNumBoundaryEdges() const;
	//! Returns the number of uknowns in the problem
	std::size_t getNumUnk() const;
	//! Returns the basis function id associated with point node
	std::size_t getNodeBF(Point2D* node);
	//! Returns element iElem
	Element<Point2D>* getElement(std::size_t iElem) const;
	//! Returns edge iEdge
	ElemEdge<Point2D>* getEdge(std::size_t iEdge) const;
	//! Returns boundary edge iEdge
	ElemEdge<Point2D>* getBoundaryEdge(std::size_t iEdge) const;
	//! Returns node iNode
	Point2D* getNode(std::size_t iNode) const;
	//@}

	//! Outputs the boundary edges to file stream outFile, using a Matlab-like format
	void printBoundaries(std::ofstream& outFile) const;
	//! Returns the total area of the mesh
	double computeArea() const;
	//! Returns the integral of material property propID over the entire mesh
	double computeMaterialPropertyTotal(unsigned propID) const;
	//! Fills inPtVec with all nodes 
	void makeNodeVec(std::vector<Point2D*>& inPtVec) const;
	//! Fills edgeVec will all edges connected to tstPt
	void getNeighboringEdges(Point2D* tstPt, std::vector<ElemEdge<Point2D>*>& edgeVec);
private:
	std::size_t numElements, numEdges, numNodes;
	std::size_t numInternalEdges, numBoundaryEdges, numUnk;
	std::vector<std::unique_ptr<Point2D>> nodeVec;
	std::vector<std::unique_ptr<ElemEdge<Point2D>>> edgeVec;
	std::vector<std::unique_ptr<Element<Point2D>>> elemVec;
	std::map<Point2D*, std::size_t> basisNodeMap;

	void finishSetup(bool& err, const std::vector<Point2D*>& inNodeVec, const std::vector<Element<Point2D>*>& inElemVec);
	bool setConnectivity();
	void putBoundaryEdgesLast();
	void ensureConsistentNormals();
	void makeUpwardNormals();
	void setStandardEdgeOrder();
	void setNodeBFMap();
	void setPatch2BasisMap();
	void readPatch(std::ifstream& geoFile, unsigned int patNum);
	bool cleanSurfaceDefinition(double& area);
	void outputInfo(double area) const;
	std::size_t findNode(const Point2D* theNode) const;
	Point2D* findBoundaryPoint(Point2D pt);
	void reorderRCM();
	int getNodeDegree(Point2D* tstPt);
	void insertSortDegree(Point2D* inPt, std::list<Point2D*>& levelSet, 
						  std::list<int>& levelSetDegrees);
	std::unique_ptr<Element<Point2D>> copyElement(const Element<Point2D>* const inElem, const std::vector<std::size_t>& ptIDVec) const;
};

inline
void Mesh2D::makeNodeVec(std::vector<Point2D*>& inPtVec) const
{
	inPtVec.clear();
	for(std::size_t k = 0; k < numNodes; k++)
		inPtVec.push_back(nodeVec[k].get());
}

inline 
std::size_t Mesh2D::getNumEdges() const 
{
	return numEdges;
}

inline 
std::size_t Mesh2D::getNumNodes() const 
{
	return numNodes;
}

inline 
std::size_t Mesh2D::getNumInternalEdges() const 
{
	return numInternalEdges;
}

inline 
std::size_t Mesh2D::getNumBoundaryEdges() const 
{
	return numBoundaryEdges;
}

inline 
std::size_t Mesh2D::getNumUnk() const 
{
	return numUnk;
}

inline 
Element<Point2D>* Mesh2D::getElement(std::size_t iElem) const 
{
	return elemVec[iElem].get();
}

inline 
ElemEdge<Point2D>* Mesh2D::getEdge(std::size_t iEdge) const 
{
	return edgeVec[iEdge].get();
}

inline
ElemEdge<Point2D>* Mesh2D::getBoundaryEdge(std::size_t iEdge) const
{
    return edgeVec[iEdge + numInternalEdges].get();
}

inline
Point2D* Mesh2D::getNode(std::size_t iNode) const
{
	return nodeVec[iNode].get();
}

inline
std::size_t Mesh2D::getNodeBF(Point2D* node)
{
	return basisNodeMap[node];
}
#endif
