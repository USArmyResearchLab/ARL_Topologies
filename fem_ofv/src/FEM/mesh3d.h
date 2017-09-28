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

#ifndef MESH3D_H
#define MESH3D_H
#include "UTIL/topologiesdefs.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <stack>
#include <queue>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>
#include <Eigen/Dense>
#include "point3d.h"
#include "element.h"
#include "elemedge.h"
#include "femmesh.h"
#include "cell.h"
#include "UTIL/genericmaterial.h"

namespace Topologies{
class TOMesh;
}

//! A volume object that completely describes one connected volume.  
/*!
Mesh3D contains all cells (which are made up of patches and nodes) and connectivity information. It can also access element matrix construction functions that are defined in the element types.  This class basically contains all necessary geometric, topological, and numerical information to describe a finite element problem.
*/
class Mesh3D : public FEMMesh
{
public: 
	//! @name Constructors and destructor
	//@{
	//! Constructor which takes a TOMesh and copies the nodes and connectivity from it
	Mesh3D(const Topologies::TOMesh* const inMesh, const Topologies::GenericMaterial& baseMat);
	//! Constructor which takes in a vector of Point3D pointers and a vector of Cell objects. err is returned with any error conditions
	/*! This constructor generates a 3D mesh from a set of nodes (inNodeVec) and elements (inElemVec).  These objects will be copied into for internal use and the mesh will be processed to generate and store topological information such as cell edges and faces.  err indicates whether an error has occured on return during the mesh processing.
	*/
	Mesh3D(bool& err, const std::vector<Point3D*>& inNodeVec, const std::vector<Cell*>& inElemVec);
	Mesh3D(const Mesh3D&);
	Mesh3D(Mesh3D&&);
	Mesh3D& operator=(Mesh3D);
	void swap(Mesh3D&);
	virtual ~Mesh3D();
	//@}

	//! @name FEM element matrix functions
	//@{
	//! Returns the linear elastic element matrix
	virtual Eigen::MatrixXd getElementMatrix(std::size_t k) const {return getCell(k)->getElemMat();}
	//! Returns the linear elastic element matrix
  virtual Eigen::MatrixXd getElementMatrix(std::size_t k, const ElasticProblemType theEPT) const {return getCell(k)->getElemMat();}
	//! Returns the element matrix for the Laplacian operator for element k and material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(std::size_t k, double matVal) const {return getCell(k)->getLaplacianElemMat(matVal);}
	//! Returns the element matrix for the Laplacian operator for element k and anisotropic material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(std::size_t k, const Eigen::MatrixXd& matVal) const {return getCell(k)->getLaplacianElemMat(matVal);}
	//! Returns the stiffness matrix associated with isotropic thermal expansion
	virtual Eigen::MatrixXd getThermalExpansionMat(std::size_t k, double alphaTE, const ElasticProblemType inEPT = eptPlaneStress) const {return getCell(k)->getThermalExpansionMat(alphaTE);}
	//! Returns a vector containing the basis function values at the given local are coordinate (xi, eta, zeta)
	virtual Eigen::VectorXd getInterpolationVec(std::size_t k, double xi, double eta, double zeta) const {return getCell(k)->getInterpolationVec(xi, eta, zeta);}
	//! Returns a matrix giving the gradient of all basis functions on Cell k at the centroid of the element
	virtual Eigen::MatrixXd getGradBF(std::size_t k) const {return getCell(k)->getGradBF();}
	//! Returns a matrix giving the Laplacian of the basis functions, not implemented in 3d
	virtual Eigen::MatrixXd getLaplacianBF1(std::size_t k) const {return getCell(k)->getLaplacianBF1();}
	//! Returns a matrix giving the Laplacian of the basis functions, not implemented in 3d
	virtual Eigen::MatrixXd getLaplacianBF2(std::size_t k) const {return getCell(k)->getLaplacianBF2();}
	//! Returns the mass matrix
	/*! The mass matrix is the inner product of the basis & testing functions over an element, integrated with the density.
	*/
	virtual Eigen::MatrixXd getMassMatrix(std::size_t k) const {return getCell(k)->getMassMat();}
	//! Returns a matrix with the integral of all basis functions over Cell k
	virtual Eigen::MatrixXd getRHSTestingMatrix(std::size_t k) const {return getCell(k)->getTFInt();}
	//@}
	//! Output the mesh to a file with name fileName.  Uses a Matlab-like format
	virtual void printMesh(const std::string& fileName) const;
	//! Returns a copy of the mesh wrapped in a unique_ptr
	virtual std::unique_ptr<FEMMesh> clone() const {return std::unique_ptr<FEMMesh>(new Mesh3D(*this));}
	//! Returns the dimension of the mesh
	virtual unsigned getDim() const {return 3;}
	//! Returns the material properties associated with Cell k
	virtual Topologies::GenericMaterial getGenericMaterial(std::size_t k) const {return getCell(k)->getGenericMaterial();}

	//!@name Data access, accesible from abstract class
	//@{
	//! Returns the number of elements (Cell objects) in the mesh
	virtual std::size_t getNumElements() const {return itsNumCells;}
	//! Returns the number of nodes (or unknowns) in the mesh
	virtual std::size_t getNumUnknowns() const {return itsNunk;}
	//! Returns the global basis function id for node klocbf on Cell kelem
	virtual std::size_t getGlobalBF(std::size_t kelem, unsigned klocbf) const {return getCell(kelem)->getGlobalBF(klocbf);}
	//! Returns the number of nodes associated with element (Cell) k
	virtual unsigned getNumElementNodes(std::size_t k) const {return getCell(k)->getNumNodes();}
	//@}
	//!@name Data access, accesible only from Mesh3D
	//! Returns the number of elements (Cell objects) in the mesh
	std::size_t getNumCells() const;
	//! Returns the number of patches (Cell faces) in the mesh
	std::size_t getNumPatches() const;
	//! Returns the number of nodes (or points) in the mesh
	std::size_t getNumNodes() const;
	//! Returns the number of non-boundary Cell faces
	std::size_t getNumInternalPatches() const;
	//! Returns the number of boundary Cell faces
	std::size_t getNumBoundaryPatches() const;
	//! Returns the number of unknowns or nodes
	std::size_t getNunk() const;
	//! Returns the basis function id offset, for problems with multiple meshes
	std::size_t getOffset() const;
	//! Sets the basis function id offset
	void setOffset(std::size_t off);
	//! Returns Cell iCell
	Cell* getCell(std::size_t iCell) const;
	//! Returns patch iPatch
	Element<Point3D>* getPatch(std::size_t iPatch) const;
	//@}
	//! Returns the volume of the mesh (sum of all element volumes)
	double computeVolume() const;
	//! Returns the surface area of the bounding volume
	double computeArea() const;
private:
	// Data Members
	std::size_t itsNumCells, itsNumPatches, itsNumEdges, itsNumNodes;
	std::size_t itsFaceNunk, itsInternalNunk;
	unsigned typeOfCell;
	std::size_t itsNumInternalPatches, itsNumBoundaryPatches, itsNunk, itsOffset;
	std::vector<std::unique_ptr<Point3D>> nodeVec;
	std::vector<std::unique_ptr<ElemEdge<Point3D>>> edgeVec;
	std::vector<std::unique_ptr<Element<Point3D>>> patchVec;
	std::vector<std::unique_ptr<Cell>> cellVec;
	std::map<Point3D*, std::size_t> basisNodeMap;

	// Private Methods
	void finishSetup(bool& err, const std::vector<Point3D*>& inNodeVec, const std::vector<Cell*>& inElemVec);
	std::unique_ptr<Cell> copyCell(const Cell* const inCell, const std::vector<std::size_t>& ptIDVec) const;
	bool setPatchConnectivity();
	bool setEdgeConnectivity();
	void putBoundaryPatchesLast();
	void ensurePositiveJacobian();
	void setStandardPatchOrder();
	void setPatchBFMap();
	void setCellBFMap();
	void setNodeBFMap();
	void setCell2BasisMap();
	void cleanVolumeDefinition(double& area, double& volume);
	void outputInfo(double area, double volume) const;
	std::size_t findCell(const Cell* theCell) const;
	std::size_t findNode(const Point3D* theNode) const;
};

inline std::size_t Mesh3D::getNumCells() const
{
	return itsNumCells;
}
inline std::size_t Mesh3D::getNumPatches() const
{
	return itsNumPatches;
}
inline std::size_t Mesh3D::getNumNodes() const 
{
	return itsNumNodes;
}
inline std::size_t Mesh3D::getNumInternalPatches() const 
{
	return itsNumInternalPatches;
}
inline std::size_t Mesh3D::getNumBoundaryPatches() const 
{
	return itsNumBoundaryPatches;
}
inline std::size_t Mesh3D::getNunk() const 
{
	return itsNunk;
}
inline std::size_t Mesh3D::getOffset() const 
{
	return itsOffset;
}
inline Cell* Mesh3D::getCell(std::size_t iCell) const 
{
	return cellVec[iCell].get();
}
inline Element<Point3D>* Mesh3D::getPatch(std::size_t iPatch) const 
{
	return patchVec[iPatch].get();
}

#endif
