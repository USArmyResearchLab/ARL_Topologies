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

#ifndef CELL_H
#define CELL_H

#include <cmath>
#include <iostream>
#include <cassert>
#include <bitset>
#include <vector>
#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "element.h"
#include "lintri.h"
#include "UTIL/genericmaterial.h"
#include "point3d.h"

enum CellType
{
	ctLinTetra = 1,
	ctQuadTetra,
	ctTriLinHex,
	ctSphTetra = 100
};

//! Cell is an abstract base class with derived classes that represent different types of canonical cell elements in 3d
/*! A union of (connected) cells is what makes up a volume.   Currently tetrahedral and hexahedral elements are implemented, 
 * each of which have their own derived classes.  Some functions that compute geometric quantities are based on the paper 
 * Graglia, R.D., et. al. "Higher Order Interpolatory Vector Bases for CEM,"  IEEE TAP vol.45 pp.329-342.  (GWP from now on)
 */
class Cell
{
public:
	//! Constructor that takes a GenericMaterial
	Cell(CellType inCellType, const Topologies::GenericMaterial& inMat);
	virtual ~Cell() {};

	//! @name Data access
	//@{
	//! Returns the number of nodes associated with this Cell
	virtual unsigned short int getNumNodes() const = 0;
	//! Returns the number of patches (faces) associated with this Cell
	virtual unsigned getNumPatches() const = 0;
	//! Returns the patch (face) with id patchNum
	virtual Element<Point3D>* getPatch(unsigned short int patchNum) const = 0;
	//! Returns the node (point) with id nodeNum
	virtual Point3D* getNode(unsigned short int nodeNum) const = 0;
	//! Returns the global basis function id with local id BFNum
	std::size_t getGlobalBF(unsigned short int BFNum) const;
	//! Returns the local basis function id with global id BFNum
	unsigned short int getLocalBF(unsigned BFNum) const;
	//! Returns the cell type
	CellType getCellType() const {return itsCellType;}
	//! Return whether or not this Cell contains node inNode
	virtual bool containsNode(Point3D* inNode) const = 0;
	//! Returns the GenericMaterial associated with this Cell
	Topologies::GenericMaterial getGenericMaterial() const {return itsMaterial;}
	//@}
	//! @name Element matrix computation functions
	//@{
	//! Computes and returns the linear elastic element matrix
	virtual Eigen::MatrixXd getElemMat() const = 0;
	//! Computes and returns the Laplacian element matrix for given material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(double matVal) const = 0;
	//! Computes and returns the Laplacian element matrix for given anisotropic material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(const Eigen::MatrixXd& matVal) const = 0;
	//! Computes and returns the mass matrix (integral of basis * test functions)
	virtual Eigen::MatrixXd getMassMat() const = 0;
	//! Computes and returns the stiffness matrix associated with isotropic thermal expansion
	virtual Eigen::MatrixXd getThermalExpansionMat(double alphaTE) const = 0;
	//! Computes and returns the integral of the testing functions
	virtual Eigen::MatrixXd getTFInt() const = 0;
	//! Returns a vector containing the basis function values at the given local area coordinate (\xi, \eta, \zeta)
	virtual Eigen::VectorXd getInterpolationVec(double xi, double eta, double zeta) const = 0;
	//! Returns a matrix containing the gradient of the basis functions at element k, computed at the centroid
	virtual Eigen::MatrixXd getGradBF() const = 0;
	//! Computes and returns the Laplacian of the basis function: div(grad(phi))
	virtual Eigen::MatrixXd getLaplacianBF1() const = 0;
	//! Computes and returns the Laplacian of the basis function: grad(grad(phi))
	virtual Eigen::MatrixXd getLaplacianBF2() const = 0;
	//@}
	//! @name Geometric functions
	//@{
	//! Computes edge vectors for a Cell (see GWP) at the point (\xi_1, \xi_2, \xi_3)
	virtual void justL(double xi1, double xi2, double xi3, Point3D& l1, Point3D& l2, Point3D& l3) const = 0;
	//! Compute the edge vectors along with the Jacobian of the Cell at the point (\xi_1, \xi_2, \xi_3)
	virtual void jacL(double xi1, double xi2, double xi3, double& jac,Point3D& l1, Point3D& l2, Point3D& l3) const = 0;
	//! Returns the position in real space at the point (\xi_1, \xi_2, \xi_3)
	virtual Point3D pos(double xi1, double xi2, double xi3) const = 0; 
	//! Returns the volume of the Cell
	virtual double volumeIntegral() const = 0;
	//@}
	//! @name Topological functions
	//@{
	//! Returns whether or not the cell is on the boundary of its containing volume
	/*! A cell is on the boundary of an element if one of its patches (faces) is connected to only 1 cell */
	virtual bool isBoundaryCell() const = 0;
	//! Swaps vectors that define the area coordinates.
	/*! This is to ensure positivity of the Jacobian */
	virtual void switchXi2Xi3() = 0;
	//! Return a patch (face) of the Cell with id patchNum
	/*! The patch is defined as an Element in 3-space, and is returned in a unique_ptr */
	virtual std::unique_ptr<Element<Point3D>> createPatch(unsigned short int patchNum) const = 0;
	//! Associates Element thePatch with this Cell
	virtual void addPatch(Element<Point3D>* thePatch) = 0;
	//! Reorders the patches associated with this Cell so that they are correctly associated with given area coordinates
	virtual void setStandardPatchOrder() = 0;
	//! Adds basis function ids for a given patch theCPatch.
	/*! Patch basis functions have interpolation points on the face of the given patch */
	virtual void addPatchBases(Element<Point3D>* theCPatch, Point3D interpPos[], std::size_t startPoint, int signum) = 0;
	//! Adds basis function ids to the Cell.
	/*! Cell basis functions are completely internal to a patch (do not exist on faces, edges, or nodes)*/
	virtual void addCellBases(std::size_t& lastAdded) = 0;
	//! Adds nodal basis function ids
	/*! Nodal basis functions are only associated with Cell nodes*/
	virtual void addNodeBases(const Point3D* theNode, std::size_t bfGID) = 0;
	//@}
protected:
	std::vector<std::size_t> bfVec;
	CellType itsCellType;
	Topologies::GenericMaterial itsMaterial;
};

inline 
std::size_t Cell::getGlobalBF(unsigned short int BFNum) const 
{
	assert(BFNum < bfVec.size());
	return bfVec[BFNum];
}

#endif
