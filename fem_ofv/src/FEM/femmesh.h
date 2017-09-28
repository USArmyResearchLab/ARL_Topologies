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

#ifndef FEMMESH
#define FEMMESH

#include <string>
#include <memory>
#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "UTIL/genericmaterial.h"
#include "globaldefs.h"

//! Abstract base class for finite element meshes.  
/*!
This class basically contains all necessary geometric, topological, and numerical information to describe a finite element problem.  It is implemented in Mesh2D and Mesh3D.
*/
class FEMMesh
{
public:
	FEMMesh(){};
	virtual ~FEMMesh(){};
	//! @name Data accessors
	//@{
	//! Returns the number of elements (Cell objects) in the mesh
	virtual std::size_t getNumElements() const = 0;
	//! Returns the number of nodes (or points)
	virtual std::size_t getNumUnknowns() const = 0;
	//! Returns the number of nodes associated with element k
	virtual unsigned getNumElementNodes(std::size_t k) const = 0;
	//! Returns the global basis function id associated with element kelem for local basis function klocbf
	virtual std::size_t getGlobalBF(std::size_t kelem, unsigned klocbf) const = 0;
	//@}
	//! @name Finite element matrix compuation functions
	//@{
	//! Returns the plane stress (for 2d), linear elastic element matrix for element k
	virtual Eigen::MatrixXd getElementMatrix(std::size_t k) const = 0;
	//! Returns either the plane strain or plane stress, linear elastic element matrix for element k
	virtual Eigen::MatrixXd getElementMatrix(std::size_t k, const ElasticProblemType theEPT) const = 0;
	//! Returns the Laplacian element matrix for given material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(std::size_t k, double matVal) const = 0;
	//! Returns the Laplacian element matrix for given anisotropic material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(std::size_t k, const Eigen::MatrixXd& matVal) const = 0;
	//! Returns the stiffness matrix associated with isotropic thermal expansion
	virtual Eigen::MatrixXd getThermalExpansionMat(std::size_t k, double alphaTE, const ElasticProblemType inEPT = eptPlaneStress) const = 0;
	//! Returns a vector containing the basis function values at the given local are coordinate (xi, eta, zeta)
	virtual Eigen::VectorXd getInterpolationVec(std::size_t k, double xi, double eta, double zeta = 0.) const = 0;
	//! Returns a matrix containing the gradient of the basis functions at element k, computed at the centroid
	virtual Eigen::MatrixXd getGradBF(std::size_t k) const = 0;
	//! Returns a matrix giving the Laplacian of the basis functions, not implemented in 3d
	virtual Eigen::MatrixXd getLaplacianBF1(std::size_t k) const = 0;
	//! Returns a matrix giving the Laplacian of the basis functions, not implemented in 3d
	virtual Eigen::MatrixXd getLaplacianBF2(std::size_t k) const = 0;
	//! Returns the mass matrix
	/*! The mass matrix is the inner product of the basis & testing functions over an element, integrated with the density.
	*/
	virtual Eigen::MatrixXd getMassMatrix(std::size_t k) const = 0;
	//! Returns the integral of the testing functions over element k
	virtual Eigen::MatrixXd getRHSTestingMatrix(std::size_t k) const = 0;
	//@}
	//! Outputs this mesh to the file named in fileName.  Uses a Matlab-like format
	virtual void printMesh(const std::string& fileName) const = 0;
	//! Returns the dimension of the mesh (2 or 3)
	virtual unsigned getDim() const = 0;
	//! Returns the material properties of element k
	virtual Topologies::GenericMaterial getGenericMaterial(std::size_t k) const = 0;
	//! Returns a copy of this FEMMesh, wrapped in a unique_ptr
	virtual std::unique_ptr<FEMMesh> clone() const = 0;
private:
};

#endif
