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

#include "topopt.h"
#include "mpihandler.h"
#include "topoptobjfun.h"

namespace Topologies{
std::pair<double, bool> TopOpt::evaluateSingleObjective(const TopOptRep* pTOR, EvalFunction inef) const
{
	std::pair<std::vector<double>, bool> res = evaluateMultiObjective(pTOR, inef);
	if(res.first.empty())
		return std::make_pair(1e6, false);
	return std::make_pair(res.first[0], res.second);
}

std::pair<std::vector<double>, bool> TopOpt::evaluateMultiObjective(const TopOptRep* pTOR, EvalFunction inef) const
{
	std::pair<std::vector<double>, bool> res;
	if(pMPIH)
		res = pMPIH->rootEvaluateTOR(pTOR, inef);
	else
	{
		if(inef == efF)
			pObjFun->f(*pTOR, res);
		else // efC
			pObjFun->c(*pTOR, res);
	}
	return res;
}

std::pair<std::vector<double>, bool> TopOpt::evaluateGradient(const TopOptRep* pTOR, EvalFunction inef) const
{
	std::pair<std::vector<double>, bool> resG;
	if(pMPIH)
		resG = pMPIH->rootEvaluateTOR(pTOR, inef);
	else
	{
		if(inef == efG)
			pObjFun->g(*pTOR, resG);
		else
			pObjFun->gc(*pTOR, resG);
	}
	return resG;
}

std::tuple<double, std::vector<double>, bool> TopOpt::evaluateSingleObjectiveAndGradient(const TopOptRep* pTOR,
																																												EvalFunction inef) const
{
	std::pair<std::vector<double>, bool> fRes, gRes;
	if(pMPIH)
	{
		fRes = pMPIH->rootEvaluateTOR(pTOR, inef);
		// Remove gradient
		assert(fRes.size() > pTOR->getDataSize());
		gRes.first.insert(gRes.first.end(), fRes.first.end() - pTOR->getDataSize(), fRes.first.end());
		fRes.first.erase(fRes.first.end() - pTOR->getDataSize(), fRes.first.end());
	}
	else
	{
		if(inef == efFandG)
			pObjFun->fAndG(*pTOR, fRes, gRes);
		else
			std::cout << "Warning: Computing C and its gradient not yet supported, use evaluateSingleObjective" << std::endl;
	}
	double fval = 1e6;
	if(!fRes.first.empty())
		fval = fRes.first[0];
	return std::make_tuple(fval, gRes.first, fRes.second && gRes.second);
}

void TopOpt::filterGradient(const std::vector<double>& x, std::vector<double>& locGrad,
															TopOptRep& workTOR, double filterSize, double minDensity) const
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

