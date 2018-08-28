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

#ifndef TOFEMOBJFUN_H
#define TOFEMOBJFUN_H

#include <string>
#include <memory>
#include "OBJFUN/topoptobjfun.h"
#include "REP/cgal_types.h"
#include "femproblem.h"
#include "IO/inputloader.h"
#include "boundarycondition.h"
#include "loadcondition.h"

namespace Topologies{class TopOptRep;}

//! Interface between finite element solver and TopOptObjFun base class
/*!
 * This class is derived from TopOptObjFun and interfaces between it and FEMProblem.  FEMProblem is a static, linear elastic solver with various element types.  This class will construct an FEMProblem object, solve a given problem, and return the compliance which is the dot product of the displacement and load vectors.
 */
class TOFEMObjFun : public Topologies::TopOptObjFun
{
public:
	//! Constructor which takes the name of an input file.
	TOFEMObjFun(const std::string& inpFileName);
	virtual ~TOFEMObjFun();

	//! Print a result file for a given TopOptRep, not implemented
	virtual void printResult(const Topologies::TopOptRep& inTOR, const std::string& fileName) const {} // Not yet implemented
	//! The parentheses operator takes a TopOptRep (a topology), runs FEMProblem, and returns the compliance.
	virtual std::pair<double, bool> operator()(const Topologies::TopOptRep& inTOR) const;
	//! Function f evaluates the objective function (same as operator()).
	virtual void f(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
	//! Function c evaluates the constraints on a given TopOptRep, here the constraint is volume fraction
	virtual void c(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
	//! Function g computes the gradient of the objective function
  virtual void g(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
	//! Computes the gradient of the constraints
	virtual void gc(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
	//! Computes the gradient and the objective function value
	virtual void fAndG(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& fRes,
											std::pair<std::vector<double>, bool>& gRes) const; 
#ifdef USE_MPI
	//! This overloaded parentheses operator takes an MPI communicator for parallelized objective functions, though it is not implemented here.
	virtual std::pair<double, bool> operator()(const Topologies::TopOptRep& inTOR, MPI::Comm& communicator) const;
	//! Overloaded f with MPI communicator
	virtual void f(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
	//! Overloaded c with MPI communicator
	virtual void c(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
	//! Function g computes the gradient of the objective function
  virtual void g(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
	//! Computes the gradient of the constraints
  virtual void gc(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
	//! Computes the gradient and the objective function
	virtual void fAndG(const Topologies::TopOptRep& inTOR, std::pair<std::vector<double>, bool>& fRes,
											std::pair<std::vector<double>, bool>& gRes, MPI::Comm& communicator) const;
#endif
private:
	std::vector<ExoBC> generateExoBCVec(const Topologies::TOMesh* const inMesh, std::size_t kLoad) const;
	std::unique_ptr<FEMProblem> setupAndSolveFEM(const Topologies::TopOptRep& inTOR, std::size_t kLoad) const;
	std::unique_ptr<Topologies::TOMesh> getTOMesh(const Topologies::TopOptRep& inTOR) const;
	void g(const Topologies::TopOptRep& inTOR, const std::vector<std::map<std::size_t, double>>& dF, 
		std::pair<std::vector<double>, bool>& outRes) const;
	unsigned dim;
	Topologies::GenericMaterial baseMat;
	double maxDisplacement, volfrac;
	std::vector<BoundaryCondition> bcVec;
	std::vector<std::vector<LoadCondition<double>>> lcVV;
};

// Factory load/destroy functions to interface with optimization code
extern "C" Topologies::TopOptObjFun* create(const std::string& inpFileName) 
{
	return new TOFEMObjFun(inpFileName);
}

extern "C" void destroy(Topologies::TopOptObjFun* p) 
{
	delete p;
}
// Functions defined in tofor:
// typedef TopOptObjFun* create_toof(const std::string& inpFileName);
// typedef void destroy_toof(TopOptObjFun*);

//! XML input parser for the FEM problem
class FEMInputLoader : public Topologies::InputLoader::InputParser
{
public:
	FEMInputLoader() : dim(2), rho(1.), E(1.), nu(0.25), maxDisplacement(1e16), volfrac(0.5), inputMFF(Topologies::mffUnknown) {}
	virtual ~FEMInputLoader() {}
	virtual void parse(const pugi::xml_document& xmldoc);
	virtual void parse(const pugi::xml_node& rootNode);
	
	//! Returns Young's modulus
	double getYoungsMod() const {return E;}
	//! Returns Poisson's ratio
	double getPoissonsRatio() const {return nu;}
	//! Returns the mesh dimensions (2 or 3)
	unsigned getDim() const {return dim;}
	//! Returns the maximum allowed displacement for a result
	double getMaxDisplacement() const {return maxDisplacement;}
	//! Returns the density
	double getDensity() const {return rho;}
	//! Returns the volume fraction
	double getVolFrac() const {return volfrac;}
	//! Returns the Exodus mesh file name
	const std::string& getMeshFileName() const {return meshFilename;}
	//!Returns the mesh file format 
	Topologies::MeshFileFormat getMeshFileFormat() const {return inputMFF;};
	//! Returns the vector of boundary condition definitions
	const std::vector<BoundaryCondition>& getBCVec() const {return bcVec;} 
	//! Returns the vector of load condition definitions
  const std::vector<std::vector<LoadCondition<double>>>& getLCVec() const {return lcVV;}
private:
	static BCType parseBCType(const std::string& inBCStr);
	static CoordinateSystem::Type parseCSType(const std::string& inCSTStr);
	void setFileFormat(const pugi::xml_node& rootNode);
	void parseExoLC(const pugi::xml_node& rootNode);
	void parseExoBC(const pugi::xml_node& rootNode);
	void parseGeoBC(const pugi::xml_node& rootNode);
	void parseGeoLC(const pugi::xml_node& rootNode);
	CoordinateSystem::Type readCoordinateSystemAttribute(const pugi::xml_node& rootNode) const;

	unsigned dim;
	double rho, E, nu, maxDisplacement;
	double volfrac;
	Topologies::MeshFileFormat inputMFF;
	std::string meshFilename;
	std::vector<BoundaryCondition> bcVec;
	std::vector<std::vector<LoadCondition<double>>> lcVV;
};

#endif

