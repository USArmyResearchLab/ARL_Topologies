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

#include <sstream>
#include "topoptoc.h"
#include "topoptobjfun.h"
#include "topoptrep.h"
#include "inputloaderopt.h"
#include "outputwriter.h"
#include "geometrytranslation.h"
#include "postprocess.h"
#include "mpihandler.h"
#include "tomesh.h"

namespace Topologies{
TopOptOC::TopOptOC(const InputLoader::TOOGeneric& inputData, TopOptObjFun* inpObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inMPIH):
	TopOpt(inpObjFun, inOutVec, inMPIH),
	filterSize(inputData.getFilterSize()), 
	stopTol(inputData.getStopTol()),
	maxIters(inputData.getMaxIters()),
	minDensity(0.001)
{
}

TopOptOC::~TopOptOC()
{
}

std::unique_ptr<TopOptRep> TopOptOC::optimize(const TopOptRep& initialGuess)
{
	double curTol = stopTol + 1.;
	unsigned curIter = 0;
	std::vector<double> x;
	initialGuess.getRealRep(x);
	std::vector<std::size_t> sizes;
	initialGuess.getDataSize(sizes);
	std::cout << "Starting OC, max iterations: " << maxIters << ", stop tolerance: " << stopTol << std::endl;
	std::cout << "  Norm of initial vector: " << HelperNS::norm(x) << std::endl;
	std::unique_ptr<TopOptRep> result = initialGuess.clone();
	while(curTol > stopTol && curIter < maxIters)
	{
		++curIter;
		result->setRealRep(x.begin(), x.end());
		std::tuple<double, std::vector<double>, bool> fAndG = evaluateSingleObjectiveAndGradient(result.get(), efFandG);
//		std::pair<double, bool> curf = evaluateSingleObjective(result.get(), efF);
		std::vector<double> g = std::move(std::get<1>(fAndG));
//		computeGradient(x, g, *result);
		filterGradient(x, g, *result, filterSize, minDensity);
		curTol = ocUpdate(x, g, *result);
		std::cout << "Done iteration #" << curIter << ", current ofv: " << std::get<0>(fAndG) 
							<< ", Norm of gradient: " << HelperNS::norm(g) 
							<< ", change: " << curTol << std::endl;
		std::cout << "  volume fraction: " << result->computeVolumeFraction() << std::endl;
		result->setRealRep(x.begin(), x.end());
		handleOutput(result.get());
	}
	if(curIter >= maxIters)
		std::cout << "Done: Max iterations exceeded, current iteration: " << curIter << std::endl;
	if(curTol <= stopTol)
		std::cout << "Done: Tolerance converged, current tolerance: " << curTol << std::endl;
	result->setRealRep(x.begin(), x.end());
	return result;
}

void TopOptOC::computeGradient(const std::vector<double>& x, std::vector<double>& g, TopOptRep& workTOR)
{
	workTOR.setRealRep(x.begin(), x.end());
	std::pair<std::vector<double>, bool> resG = evaluateGradient(&workTOR);
	// Apply gradient filter
	filterGradient(x, resG.first, workTOR, filterSize, minDensity);
	// Copy to g
	g = std::vector<double>(resG.first);
}

double TopOptOC::ocUpdate(std::vector<double>& x, const std::vector<double>& g, TopOptRep& toRep)
{
	double l1 = 0., l2 = 100000., move = 0.05;
	std::vector<double> xnew(x.size());
	while((l2 - l1) > 1e-4)
	{
		// Update xnew
		double lmid = 0.5*(l1 + l2);
		for(std::size_t k = 0; k < xnew.size(); ++k)
		{
			double updateVal = x[k]*sqrt(fabs(g[k])/lmid);
			xnew[k] = MAX(0., MAX(x[k] - move, MIN(1., MIN(x[k] + move, updateVal))));
		}
		// Recompute constraint
		toRep.setRealRep(xnew.begin(), xnew.end());
		std::pair<double, bool> curc = evaluateSingleObjective(&toRep, efC);
		double c = curc.first;
		if(!curc.second)
		{
			std::cout << "Error: Constraint computation in ocUpdate returned bad value, aborting" << std::endl;
			abort();
		}
		if(c > 0.)
			l1 = lmid;
		else
			l2 = lmid;
	}
	// Compute max change in x vals
	double maxChange = 0.;
	for(std::size_t k = 0; k < xnew.size(); ++k)
	{
		if(fabs(x[k] - xnew[k]) > maxChange)
			maxChange = fabs(x[k] - xnew[k]);
	}
	x = xnew;
	return maxChange;
}

}

