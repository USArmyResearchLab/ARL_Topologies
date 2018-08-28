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

#ifndef TOPOPTOBJFUN_H
#define TOPOPTOBJFUN_H

#include <string>
#include <vector>
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "helper.h"
#include "topologiesdefs.h"

namespace Topologies{
class TopOptRep;
class MPIHandler;

//! Abstract base class for an objective function to be used by TopOpt optimizers
/*! Implementations of this class define the problem that will be optimized.  There are no 
*  implementations in Topologies, as these are problem specific and must be implemented
*  in a shared library.  However, the Topologies repository has a basic 2d and 3d, static,
*  linear elastic solver for implementing classic topology optimization problems.
*  Implementations can be either serial or parallel. Parallel implementations must use
*  the MPI communicator provided, and not MPI::COMM_WORLD.
*  This class does implement a parallelized finite difference method for computing gradients
*  of both the objective function and constraints.  This may be overrided in derived classes
*  if more efficient or accurate methods are available.
*/
class TopOptObjFun
{
public:
	TopOptObjFun() { setEps(); }
	virtual ~TopOptObjFun() {}
	TopOptObjFun(const TopOptObjFun&) { setEps(); }
	TopOptObjFun& operator=(const TopOptObjFun&) {setEps(); return *this;}
	TopOptObjFun(TopOptObjFun &&) { setEps(); }
	TopOptObjFun& operator=(TopOptObjFun &&) {setEps(); return *this;}

	//! Computes and returns the objective function value and a flag to indicate failure of the objective function
	virtual std::pair<double, bool> operator()(const TopOptRep& inTOR) const = 0;
	//! Computes the objective function value and stores it in outRes.  Can be multiple objectives
	virtual void f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const = 0;
	//! Computes the gradient of the objective function with respect to all design variables in the TopOptRep (in serial)
	/*! This function computes the gradient of the objective function.  Note that the default implementation is a 
	 *  finite difference method and this should be overrided in an implemenation if possible.
	 */ 
	virtual void g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
	//! Computes the objective function value (stored in @param fRes) and gradient (stored in @param gRes)
	/*! Evaluates both the objective function and its gradient, useful for self-adjoint problems that
	 *  efficiently compute both simultaneously.  If not overridden, will simply call f and g separately. 
	 */
	virtual void fAndG(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& fRes, 
											std::pair<std::vector<double>, bool>& gRes) const;
	//! Computes the constraint values for a given TopOptRep
	virtual void c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const = 0;
	//! Computes the gradient of the constraint function with respect to all design variables in the TopOptRep
	/*! This function computes the gradient of the constraint function and can use the MPIHandler to parallelize.
	 *  Note that the default implementation is a finite difference method and this should be overrided in
	 *  an implemenation if possible.  The finite difference computation is parallelized.
	 */ 
	virtual void gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
#ifdef USE_MPI
	//! Computes and returns the objective function value and a flag to indicate failure of the objective function
	/*! The MPI::Comm argument can be used to parallelize the objective function computation */
	virtual std::pair<double, bool> operator()(const TopOptRep& inTOR, MPI::Comm& communicator) const = 0;
	//! Computes the objective function value and stores it in outRes.  Can be multiple objectives
	/*! The MPI::Comm argument can be used to parallelize the objective function computation */
	virtual void f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const = 0;
	//! Computes the gradient of the objective function with respect to all design variables in the TopOptRep
	/*! This function computes the gradient of the objective function and can use the MPIHandler to parallelize.
	 *  Note that the default implementation is a finite difference method and this should be overrided in
	 *  an implemenation if possible.  The finite difference computation is parallelized.
	 */
	virtual void g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
	//! Computes the objective function value (stored in @param fRes) and gradient (stored in @param gRes)
	/*! Evaluates both the objective function and its gradient, useful for self-adjoint problems that
	 *  efficiently compute both simultaneously.  If not overridden, will simply call f and g separately.
	 */
	virtual void fAndG(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& fRes,
										std::pair<std::vector<double>, bool>& gRes, MPI::Comm& communicator) const;
	//! Computes the constraint values for a given TopOptRep
	/*! The MPI::Comm argument can be used to parallelize the objective function computation */
	virtual void c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const = 0;
	//! Computes the gradient of the constraint function with respect to all design variables in the TopOptRep (in serial)
	/*! This function computes the gradient of the constraint function.  Note that the default implementation is a
	 *  finite difference method and this should be overrided in an implemenation if possible.
	 */
	virtual void gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
#endif
	//! Generates output for a given TopOptRep.  This output is determined by the implementation and is usually some kind of results file
	virtual void printResult(const TopOptRep& inTOR, const std::string& fileName) const = 0;
private:
	void forwardDiff(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPIHandler* pMPIH,EvalFunction inef) const;
	void forwardDiffMPIDiff(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPIHandler* pMPIH, EvalFunction inef) const;
	void centeredDiff(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPIHandler* pMPIH, EvalFunction inef) const;
	void setEps();
	double getEps() const {return mEps;}
	void computeGradient(EvalFunction inEF, const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
	double mEps;
};

inline
void TopOptObjFun::f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	std::pair<double, bool> tmpRes = (*this)(inTOR);
	std::vector<double> ofvVec(1, tmpRes.first);
	outRes = std::make_pair(ofvVec, tmpRes.second);
}

#ifdef USE_MPI
inline
void TopOptObjFun::f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	std::pair<double, bool> tmpRes = (*this)(inTOR, communicator);
	std::vector<double> ofvVec(1, tmpRes.first);
	outRes = std::make_pair(ofvVec, tmpRes.second);
}
#endif
} // namespace

// Load/destroy functions to interface with shared library
typedef Topologies::TopOptObjFun* create_toof(const std::string& inpFileName);
typedef void destroy_toof(Topologies::TopOptObjFun*);

#endif

