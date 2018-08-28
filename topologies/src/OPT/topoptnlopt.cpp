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

#include "topoptnlopt.h"
#include "topoptobjfun.h"
#include "topoptrep.h"
#include "inputloaderopt.h"
#include "outputwriter.h"
#include "geometrytranslation.h"
#include "postprocess.h"
#include "mpihandler.h"
#include "tomesh.h"

#include <nlopt.hpp>

namespace Topologies{
namespace
{
	// Helper functions
	nlopt::algorithm convertToNLOpt(TOOType inTOOT)
	{
		if(inTOOT == tootMMA)
			return nlopt::LD_MMA;
		else if(inTOOT == tootBFGS)
			return nlopt::LD_LBFGS;
		std::cout << "Error in TopOptNLOpt, unknown optimizer type.  Defaulting to MMA" << std::endl;
		return nlopt::LD_MMA;
	}
	std::string getOptimizerName(TOOType inTOOT)
	{
		if(inTOOT == tootMMA)
			return "MMA";
		else if(inTOOT == tootBFGS)
			return "BFGS";
		std::cout << "Error in TopOptNLOptUnc, unknown optimizer type.  Defaulting to MMA" << std::endl;
		return "MMA";
	}
	bool optimizerSupportsConstraints(TOOType inTOOT)
	{
		if(inTOOT == tootMMA)
			return true;
		return false;
	}
}

TopOptNLOpt::TopOptNLOpt(TOOType inTOOT, const InputLoader::TOOGeneric& inputData, TopOptObjFun* inpObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inMPIH):
	TopOpt(inpObjFun, inOutVec, inMPIH),
	myTOOT(inTOOT),
	curiter(0),
	filterSize(inputData.getFilterSize()), 
	stopTol(inputData.getStopTol()),
	maxIters(inputData.getMaxIters()),
	constraintPenalty(inputData.getConstraintPenalty()),
	penaltyPower(inputData.getPenaltyPower()),
	useConstraints(optimizerSupportsConstraints(inTOOT)),
	minDensity(0.001)
{
}

TopOptNLOpt::~TopOptNLOpt()
{
}

std::unique_ptr<TopOptRep> TopOptNLOpt::optimize(const TopOptRep& initialGuess)
{
	double curTol = stopTol + 1.;
	unsigned curIter = 0;
	std::vector<double> vals;
	initialGuess.getRealRep(vals);
	std::vector<std::size_t> sizes;
	initialGuess.getDataSize(sizes);
	std::cout << "Starting " << getOptimizerName(myTOOT) << ", max iterations: " << maxIters << ", stop tolerance: " << stopTol << std::endl;
	workTOR = initialGuess.clone();

	// Set up NLOpt data structures
	nlopt::opt opt(convertToNLOpt(myTOOT), vals.size());
	// Set lower/upper bound constraints
	std::replace_if(vals.begin(), vals.end(), HelperNS::greaterThan1, 1.);
	std::replace_if(vals.begin(), vals.end(), HelperNS::lessThan0, 0.);
	std::vector<double> lb(vals.size(), 0.), ub(vals.size(), 1.);
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	// Set objective function
	curiter = 0;
	opt.set_min_objective(TopOptNLOpt::nloptOFWrapper, this);
	// Set constraint
	// Constraint uses a tolerance, consider making this an input parameter
	if(useConstraints)
		opt.add_inequality_constraint(TopOptNLOpt::nloptConstraintWrapper, this, 1e-4);
//	opt.add_equality_constraint(TopOptNLOpt::nloptConstraintWrapper, this, 1e-8);
	// Set stopping criterion
	opt.set_xtol_rel(stopTol);
	opt.set_maxeval(maxIters);
	// Optimize
	double minf;
	nlopt::result nlres = opt.optimize(vals, minf);
	std::cout << "NLOpt complete, minf: " << minf << std::endl;
	setTOR(vals);
	return std::move(workTOR);
}

void TopOptNLOpt::computeGradient(const std::vector<double>& x, std::vector<double>& g, double curc) const
{
	setTOR(x);
	std::pair<std::vector<double>, bool> resG = evaluateGradient(workTOR.get());
	g = std::move(resG.first);
	addGradientConstraints(g, curc);
	filterGradient(x, g, *workTOR, filterSize, minDensity);
}

void TopOptNLOpt::addGradientConstraints(std::vector<double>& resG, double curc) const
{
	if(!useConstraints)
	{
		// Compute penalty function gradient
		std::pair<std::vector<double>, bool> curgc = evaluateGradient(workTOR.get(), efGC);
		double fact = constraintPenalty*penaltyPower*pow(fabs(curc), penaltyPower - 1.);
		if(curc < 0.)
			fact = 0.;
//			fact *= -1.;
		for(std::size_t k = 0; k < resG.size(); ++k)
			resG[k] += fact*curgc.first[k];
	}
}

double TopOptNLOpt::nloptOFWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
	reinterpret_cast<TopOptNLOpt*>(data)->incrementCurIter();
	return reinterpret_cast<TopOptNLOpt*>(data)->f(x, grad);
}

double TopOptNLOpt::f(const std::vector<double>& x, std::vector<double>& grad) const
{
	setTOR(x);
	double curf;
	if(!grad.empty())
	{
		std::tuple<double, std::vector<double>, bool> res = evaluateSingleObjectiveAndGradient(workTOR.get(), efFandG);
		curf = std::get<0>(res);
		grad = std::move(std::get<1>(res));
	}
	else
	{
		std::pair<double, bool> res = evaluateSingleObjective(workTOR.get(), efF);
		curf = res.first;
	}
	double c0 = addConstraints(curf);
	if(!grad.empty())
	{
		addGradientConstraints(grad, c0);
		filterGradient(x, grad, *workTOR, filterSize, minDensity);
	}
	std::cout << curiter << ": Obj fun val: " << curf << std::endl;
	handleOutput(workTOR.get());
	return curf;
}

double TopOptNLOpt::addConstraints(double& curf) const
{
	// Add constraints if necessary
	double c0 = 0.;
	if(!useConstraints)
	{
		// Add constraint violation penalty
		std::pair<std::vector<double>, bool> curc = evaluateMultiObjective(workTOR.get(), efC);
		double addTerm = 0.;
		for(std::size_t k = 0; k < curc.first.size(); ++k)
		{
			if(curc.first[k] > 0.) // Inequality constraint
				addTerm += pow(curc.first[k], penaltyPower)*constraintPenalty;
		}
		curf += addTerm;
		if(!curc.first.empty())
			c0 = curc.first[0];
	}
	return c0;
}

double TopOptNLOpt::nloptConstraintWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
	return reinterpret_cast<TopOptNLOpt*>(data)->c(x, grad);
}

double TopOptNLOpt::c(const std::vector<double>& x, std::vector<double>& grad) const
{
	setTOR(x);
	// Compute constraints
	std::pair<std::vector<double>, bool> curc = evaluateMultiObjective(workTOR.get(), efC);
	// Compute gradient of constraints
	if(!grad.empty())
	{
		std::pair<std::vector<double>, bool> curgc = evaluateGradient(workTOR.get(), efGC);
		grad = std::move(curgc.first);
	}
	if(curc.first.size() > 0)
	{
		std::cout << "Constraint value: " << curc.first[0] << std::endl;
		return curc.first[0];
	}
	std::cout << "Warning: No constraints passed from objective function!" << std::endl;
	return 0;
}

void TopOptNLOpt::setTOR(const std::vector<double>& x) const
{
  workTOR->setRealRep(x);
}

}

