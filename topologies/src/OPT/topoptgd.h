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

#ifndef TOPOPTGD_H
#define TOPOPTGD_H

#include "topopt.h"
#include <memory>
#include <vector>

namespace Topologies{
class TopOptRep;
class TopOptObjFun;
class MPIHandler;
class OutputHandler;

//! Implementation of a basic gradient descent method 
/*! This is a basic gradient descent method and may not be useful in general.  
 *  It's main purpose is to be a simple method for debugging.
 */
class TopOptGD : public TopOpt
{
public:
	TopOptGD(const InputLoader::TOOGeneric& inputData, TopOptObjFun* inpObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inMPIH = nullptr);
	virtual ~TopOptGD();

	virtual std::unique_ptr<TopOptRep> optimize(const TopOptRep& initialGuess);

private:
	double gdUpdate(std::vector<double>& x, const std::vector<double>& g, double prevOFV, TopOptRep& workTOR);
	void computeGradient(const std::vector<double>& x, std::vector<double>& g, TopOptRep& workTOR);
	void filterGradient(const std::vector<double>& x, std::vector<double>& locGrad, TopOptRep& workTOR) const;

	double filterSize, minDensity, stopTol, gamma, maxStep;
	unsigned maxIters;
};
}
#endif

