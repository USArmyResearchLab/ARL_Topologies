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

#include "topoptgd.h"
#include "topoptobjfun.h"
#include "topoptrep.h"
#include "inputloaderopt.h"
#include "outputwriter.h"
#include "geometrytranslation.h"
#include "postprocess.h"
#include "mpihandler.h"
#include "tomesh.h"

namespace Topologies{
TopOptGD::TopOptGD(const InputLoader::TOOGeneric& inputData, TopOptObjFun* inpObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inMPIH):
	TopOpt(inpObjFun, inOutVec, inMPIH),
	filterSize(inputData.getFilterSize()), 
	stopTol(inputData.getStopTol()),
	maxIters(inputData.getMaxIters()),
	minDensity(0.001),
  gamma(inputData.getStepSize()),
	maxStep(inputData.getMaxStep())
{
}

TopOptGD::~TopOptGD()
{
}

std::unique_ptr<TopOptRep> TopOptGD::optimize(const TopOptRep& initialGuess)
{
	double curTol = stopTol + 1.;
	unsigned curIter = 0;
	std::vector<double> x;
	initialGuess.getRealRep(x);
	std::vector<std::size_t> sizes;
	initialGuess.getDataSize(sizes);
	std::cout << "Starting GD, max iterations: " << maxIters << ", stop tolerance: " << stopTol << std::endl;
	std::cout << "  Norm of initial vector: " << HelperNS::norm(x) << std::endl;
	std::unique_ptr<TopOptRep> workTOR = initialGuess.clone();
	std::unique_ptr<TopOptRep> bestRes = initialGuess.clone();
	// Initial ofv
	std::pair<double, bool> curf = evaluateSingleObjective(workTOR.get(), efF);
	double bestOFV = curf.first;
	while(curTol > stopTol && curIter < maxIters)
	{
		++curIter;
		std::vector<double> g(x.size(), 0.);
		computeGradient(x, g, *workTOR);
		curTol = gdUpdate(x, g, curf.first, *workTOR);
		workTOR->setRealRep(x);
		std::pair<double, bool> curf = evaluateSingleObjective(workTOR.get(), efF);
		if(curf.first < bestOFV)
		{
			bestRes = std::move(workTOR->clone());
			bestOFV = curf.first;
		}
		handleOutput(workTOR.get());
		std::cout << "Done iteration #" << curIter << ", current ofv: " << curf.first << ", Norm of gradient: " << HelperNS::norm(g)
              << ", change: " << curTol << std::endl;
	}
	if(curIter >= maxIters)
		std::cout << "Done: Max iterations exceeded, current iteration: " << curIter << std::endl;
	if(curTol <= stopTol)
		std::cout << "Done: Tolerance converged, current tolerance: " << curTol << std::endl;
	return bestRes;
}

void TopOptGD::computeGradient(const std::vector<double>& x, std::vector<double>& g, TopOptRep& workTOR)
{
	workTOR.setRealRep(x);
	std::pair<std::vector<double>, bool> resG = evaluateGradient(&workTOR);
	// Apply gradient filter
	filterGradient(x, resG.first, workTOR);
	// Copy to g
	g = resG.first;
}

double TopOptGD::gdUpdate(std::vector<double>& x, const std::vector<double>& g, double prevOFV, TopOptRep& workTOR)
{
	std::vector<double> xnew = HelperNS::vecMinus(x, HelperNS::vecScalarMult(gamma, g)); // x - gamma*g
	double maxChange = 0.;
	for(std::size_t k = 0; k < xnew.size(); ++k)
	{
		if(fabs(x[k] - xnew[k]) > maxStep)
		{
			double sign = xnew[k] > x[k] ? 1. : -1.;
			xnew[k] = x[k] + sign*maxStep;
		}
		if(fabs(x[k] - xnew[k]) > maxChange)
			maxChange = fabs(x[k] - xnew[k]);
	}
	x = xnew;
	return maxChange;
}

void TopOptGD::filterGradient(const std::vector<double>& x, std::vector<double>& locGrad, 
																TopOptRep& workTOR) const
{
	if(filterSize > 0.)
	{
		// Dot multiply gradient and x
		std::vector<double> workVec(locGrad.size());
		for(std::size_t k = 0; k < locGrad.size(); ++k)
		{
			if(x[k] < minDensity)
				workVec[k] = locGrad[k]*minDensity;
			else
				workVec[k] = locGrad[k]*x[k];
		}
		workTOR.filterData(workVec, filterSize);
		// Normalize
		for(std::size_t k = 0; k < locGrad.size(); ++k)
		{
			if(x[k] < minDensity)
				locGrad[k] = workVec[k]/minDensity;
			else
				locGrad[k] = workVec[k]/x[k];
		}
	}
}
}
