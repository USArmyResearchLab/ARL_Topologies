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

#ifndef TOPOPTNLOPT_H
#define TOPOPTNLOPT_H

#include "topopt.h"
#include "topoptrep.h"
#include <memory>
#include <vector>

namespace Topologies{
class TopOptObjFun;
class MPIHandler;
class OutputHandler;

//! Interface to NLOpt library
class TopOptNLOpt : public TopOpt
{
public:
  TopOptNLOpt(TOOType inTOOT, const InputLoader::TOOGeneric& inputData, TopOptObjFun* inpObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inMPIH = nullptr);
  virtual ~TopOptNLOpt();
	TopOptNLOpt(const TopOptNLOpt& copy);
	TopOptNLOpt(TopOptNLOpt&& copy) {swap(copy);}
	TopOptNLOpt& operator=(TopOptNLOpt rhs) {swap(rhs); return *this;}
	void swap(TopOptNLOpt& arg);

  virtual std::unique_ptr<TopOptRep> optimize(const TopOptRep& initialGuess);

	void incrementCurIter() {++curiter;} // For reporting inside NLOpt

private:
	// Objective function wrappers for NLOpt
  static double nloptOFWrapper(const std::vector<double> &x, std::vector<double> &grad, void* data);
  static double nloptConstraintWrapper(const std::vector<double> &x, std::vector<double> &grad, void* data);
	// Member functions
	double f(const std::vector<double> &x, std::vector<double> &grad) const;
	double c(const std::vector<double> &x, std::vector<double> &grad) const;
	// Utility functions needed for above
	double addConstraints(double& curf) const;
	void addGradientConstraints(std::vector<double>& resG, double curc) const;
  void computeGradient(const std::vector<double>& x, std::vector<double>& g, double curc) const;
	void computeGradientC(const std::vector<double>& x, std::vector<double>& g) const;
  void setTOR(const std::vector<double>& x) const;
	// Data
	TOOType myTOOT;
	unsigned curiter;
  double filterSize, minDensity, stopTol;
	double constraintPenalty, penaltyPower;
  unsigned maxIters;
	bool useConstraints;
	mutable std::unique_ptr<TopOptRep> workTOR; // This is used as temporary storage, no need to write custom copy constructors
};

inline
TopOptNLOpt::TopOptNLOpt(const TopOptNLOpt& copy) :
	TopOpt(copy),
	myTOOT(copy.myTOOT),
	curiter(copy.curiter),
	filterSize(copy.filterSize),
	minDensity(copy.minDensity),
	stopTol(copy.stopTol),
	constraintPenalty(copy.constraintPenalty),
	penaltyPower(copy.penaltyPower),
	maxIters(copy.maxIters),
	useConstraints(copy.useConstraints),
	workTOR(copy.workTOR->clone())
{
}

inline
void TopOptNLOpt::swap(TopOptNLOpt& arg)
{
	TopOpt::swap(arg);
	std::swap(myTOOT, arg.myTOOT);
	std::swap(curiter, arg.curiter);
	std::swap(filterSize, arg.filterSize);
	std::swap(minDensity, arg.minDensity);
	std::swap(stopTol, arg.stopTol);
	std::swap(constraintPenalty, arg.constraintPenalty);
	std::swap(penaltyPower, arg.penaltyPower);
	std::swap(maxIters, arg.maxIters);
	std::swap(useConstraints, arg.useConstraints);
	workTOR.swap(arg.workTOR);
}
}
#endif

