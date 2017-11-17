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

#ifndef FEMPROBLEM_H
#define FEMPROBLEM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "UTIL/topologiesdefs.h"
#include "femmesh.h"
#include "element.h"
#include "IO/exotxtmeshloader.h"
#include "coordinatesystem.h"
#include "point3d.h"

namespace Topologies{
class TOMesh;
class GenericMaterial;
}

//! A data structure to hold information for boundary conditions loaded from an Exodus II file
/*! This struct contains several variables needed to define a structural boundary condition.  Included are x, y, and z supports (displacement fixed to 0) and a load vector.  The set of nodes associated with this boundary condition is also stored.
*/
struct ExoBC
{
	ExoBC() : dim(2), isSupport(false), xsup(false), ysup(false), zsup(false), nodeSetID(0), 
													ct(CoordinateSystem::Type::cartesian) {}
	ExoBC(unsigned inDim) : dim(inDim), isSupport(false), xsup(false), ysup(false), zsup(false), 
													nodeSetID(0), ct(CoordinateSystem::Type::cartesian) {}
	unsigned dim;
	bool isSupport;
	bool xsup, ysup, zsup;
	Point3D loadVec;
	unsigned nodeSetID;
	std::vector<std::size_t> nodeIDVec;
	CoordinateSystem::Type ct;
};

//! A class that sets up and solves a static, linear elastic finite element problem
/*! This class takes a mesh and a base material, and sets up its own finite element mesh, which can compute element matrices and other values needed to solve an FEM problem.  The base material baseMat is modified by the optVal parameters in TOMesh.
*/
class FEMProblem
{
public:
	//! @name Constructors and destructor
	//@{
	//! Constructor that sets up a FEMProblem using a TOMesh
	explicit FEMProblem(const Topologies::TOMesh* const inMesh, const Topologies::GenericMaterial& baseMat);
	~FEMProblem();
	FEMProblem(const FEMProblem&);
	FEMProblem(FEMProblem&&);
	FEMProblem operator=(FEMProblem);
	void swap(FEMProblem&);
	//@}
	//! Change problem to boundary conditions specified in bcVec
	/*! This function is the main way to interact with FEMProblem.  Calling this will cause the FEMProblem objecto to recompute the FEM matrix and resolve it.
	*/
	void changeBoundaryConditionsTo(const std::vector<ExoBC>& bcVec);
	//! Returns the compliance (dot product of displacement and force) and whether or not an error occured during the solve
	std::pair<double, bool> computeCompliance();
	//! Returns a vector containing the solution
	const Eigen::VectorXd& getDisplacement() const {return *pVVec;}
	//! Returns constant access to the stiffness matrix
	const Eigen::SparseMatrix<double>& getFEMMatrix() const {return *pFEMMatrix;}
	//! Returns whether or not the last problem ran successfully
	bool validRun() const {return !invalid;}
	//! Returns the gradient of the compliance
	std::vector<double> gradCompliance(const Topologies::TOMesh* const inMesh, const std::vector<std::map<std::size_t, double>>& dTOR) const;
private:
	typedef Eigen::MatrixXd EigenDenseMat;
	typedef Eigen::VectorXd EigenVector;
	typedef Eigen::SparseMatrix<double> EigenSparseMat;
	typedef Eigen::Triplet<double> EigenT;

	void solveProblem();
	void setMatrix();
	void setVector();
	void assembleMatrix(std::vector<EigenT>& rseMat, const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const;
	void assembleMatrix2D(std::vector<EigenT>& rseMat, const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const;
	void assembleMatrix3D(std::vector<EigenT>& rseMat, const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const;
	double elementCompliance(const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const;
	double elementCompliance2D(const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const;
	double elementCompliance3D(const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const;
	bool checkForSimplex() const;
	Point3D getMeshPoint(std::size_t kn) const;

	std::unique_ptr<FEMMesh> probMesh;
	std::unique_ptr<EigenSparseMat> pFEMMatrix;
	std::unique_ptr<EigenVector> pVVec, pForce;
	std::vector<bool> fixedDOFs;
	std::vector<std::size_t> bfRemapVec;
	std::map<std::size_t,double> loadVec;
	unsigned dim;
	std::size_t numFreeDOFs;
	bool invalid;
	double itTol;
};

#endif
