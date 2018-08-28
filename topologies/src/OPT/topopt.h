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

#ifndef TOPOPT_H
#define TOPOPT_H

#include "inputloaderopt.h"
#include "outputhandler.h"
#include <vector>
#include <tuple>

namespace Topologies{
class MPIHandler;
class TopOptRep;
class TopOptObjFun;

//! Base class for optimization algorithms used for topology optimization
/*! This base class defines the public interface used for the optimization algorithms
*   used in Topologies.  
*/
class TopOpt
{
public:
	TopOpt() : pMPIH(nullptr), pObjFun(nullptr) {}
	//! Constructor that takes a TopOptObjFun, a vector of OutputHandler objects, and the MPIHandler
	TopOpt(TopOptObjFun* inObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inMPIH = nullptr) : 
	  pObjFun(inObjFun), outputVec(inOutVec), pMPIH(inMPIH){}
	virtual ~TopOpt() {}
	// Note: Shallow copies are used as this does not own the memory for its pointers
	void swap(TopOpt& arg);

	//! Main optimization function, takes an initial guess and returns a result TopOptRep
	/*! This is the main optimization routine.  The initial guess is set up in TopOptUniverse
	 *  (as defined by the input file), and the optimizer is run as specified in an input file.
	 *  The final result after termination of the algorithm is returned.
	 */
	virtual std::unique_ptr<TopOptRep> optimize(const TopOptRep& initalGuess) = 0;

protected:
	std::pair<double, bool> evaluateSingleObjective(const TopOptRep* pTOR, EvalFunction inef = efF) const;
	std::pair<std::vector<double>, bool> evaluateMultiObjective(const TopOptRep* pTOR, EvalFunction inef = efF) const; 
	std::pair<std::vector<double>, bool> evaluateGradient(const TopOptRep* pTOR, EvalFunction inef = efG) const;
	std::tuple<double, std::vector<double>, bool> evaluateSingleObjectiveAndGradient(const TopOptRep* pTOR, 
																																										EvalFunction inef = efFandG) const;
	void filterGradient(const std::vector<double>& x, std::vector<double>& locGrad,
											TopOptRep& workTOR, double filterSize, double minDensity) const;
	void handleOutput(const TopOptRep* const torToPrint) const;

	MPIHandler* pMPIH;
	TopOptObjFun* pObjFun;
	std::vector<OutputHandler*> outputVec;
};

inline
void TopOpt::swap(TopOpt& arg)
{
	std::swap(pMPIH, arg.pMPIH);
	std::swap(pObjFun, arg.pObjFun);
	outputVec.swap(arg.outputVec);
}

inline
void TopOpt::handleOutput(const TopOptRep* const torToPrint) const
{
	for(unsigned k = 0; k < outputVec.size(); ++k)
		outputVec[k]->handleOutput(torToPrint, pObjFun, false);
}
}
#endif

