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

#include "topoptga.h"
#include "helper.h"
#include "topoptrep.h"
#include "mpihandler.h"
#include "topoptobjfun.h"
#include "outputwriter.h"
#include "geneticoperators.h"
#include "selection.h"
#include "tomesh.h"
#include <algorithm>
#include <fstream>
#include <numeric>

namespace Topologies{
TopOptRealGA::TopOptRealGA(const InputLoader::TOOGA& inputData, TopOptObjFun* inObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inpMPIH):
	TopOpt(inObjFun, inOutVec, inpMPIH),
	constraintPenalty(inputData.getConstraintPenalty()),
	penaltyPower(inputData.getPenaltyPower()),
	mutationRange(inputData.getMutationRange()),
	crossRate(inputData.getCrossoverRate()),
	mutationRate(inputData.getMutationRate()),
	useDiscreteGA(false), // Consider making useDiscreteGA an input option, though needs further testing
	popSize(inputData.getPopSize()),
	numGens(inputData.getNumGens()),
	numGoals(inputData.getIsPareto() ? inputData.getNumGoals() : 1),
	ntourn(inputData.getNumTourn()),
	goalWeights(inputData.getGoalWeights()),
	sharingRadius(inputData.getSharingRadius()),
	torTemplate(nullptr),
	maxMutRad(inputData.getMutationRadius()),
	curMutRad(inputData.getMutationRadius()),
	numElite(inputData.getNumElite()),
	popAvg(0.),
	popBest(1e30),
	bestEver(1e30),
	popStdDev(0.),
	numValid(0)
{
}

TopOptRealGA::~TopOptRealGA()
{
}

void TopOptRealGA::initializePopulation(TopOptRep& initialGuess)
{
	if(popSize == 0)
		return;
	std::vector<double> realRep;
	initialGuess.getRealRep(realRep); // Copy realRep for later
	populationVec.resize(popSize);
	objFunValVec.resize(popSize);
	// Copy initial guess to chromo 0
	populationVec[0].resize(realRep.size());
	std::copy(realRep.begin(), realRep.end(), populationVec[0].begin());
	for(unsigned kpop = 1; kpop < popSize; ++kpop)
	{
		// Randomize: (may be used rather than the loop of mutations below)
		std::vector<double> tmpRealRep;
		initialGuess.randomize();
		initialGuess.getRealRep(tmpRealRep);
		std::list<double>& curChromo = populationVec[kpop];
		curChromo.resize(tmpRealRep.size());
		std::copy(tmpRealRep.begin(), tmpRealRep.end(), curChromo.begin());
		objFunValVec[kpop] = std::make_pair(std::vector<double>(), false);
	}
	// Mutate to start
//	for(unsigned k = 0; k < 100; ++k)
//		doMutation();
	bestEver = 1e30;
	initialGuess.setRealRep(realRep); // restore
}

std::unique_ptr<TopOptRep> TopOptRealGA::optimize(const TopOptRep& initialGuess)
{
	torTemplate = initialGuess.clone();
	initialGuess.getDataSize(chromoSizes);
	initializePopulation(*torTemplate);
	std::ofstream statFile("stat.txt");
	std::cout << "Starting GA: " << std::endl;
	bool done = false;
	for(unsigned kgen = 0; kgen < numGens && !done; ++kgen)
	{
		evaluatePop();
		done = noneValid();
		evaluateStatistics();
		printStatistics(statFile);
		outputResults();
		// Genetic operators
		doSelection();
		doElitism();
		doCrossover();
		// Recompute nonlocal mutation radius, linearly decreases with number of gens.
		double tmp = (double)maxMutRad - (double)maxMutRad/((double)numGens)*(double)kgen;
		curMutRad = HelperNS::round(tmp);
		std::cout << "mutation radius: " << curMutRad << std::endl;
		doMutation();
		boundsCheckChromos();
	}
	std::cout << "GA finishing, last evaluation" << std::endl;
	evaluatePop();
	evaluateStatistics();
	printStatistics(statFile);
	if(noneValid())
		std::cout << "GA optimizer is stopping because there are no valid chromosomes" << std::endl;
	else
		std::cout << "GA optimizer complete, best ofv: " << bestEver << std::endl;
	std::unique_ptr<TopOptRep> result = initialGuess.clone();
	if(!bestChromoVec.empty())
		setTOR(bestChromoVec[0], *result);
	statFile.close();
	return result;
}

void TopOptRealGA::evaluatePop()
{
	if(useDiscreteGA)
		makeDiscreteChromos();
	if(pMPIH)
	{
		// Set up vector of TORs for evaluation
		std::vector<std::unique_ptr<TopOptRep> > torUPVec;
		std::vector<TopOptRep*> torVec(popSize);
		for(unsigned k = 0; k < popSize; ++k)
		{
			setTOR(populationVec[k], *torTemplate);
			torUPVec.push_back(torTemplate->clone());
			torVec[k] = torUPVec[k].get();
		}
		pMPIH->rootBatchEvaluate(torVec, objFunValVec);
		addConstraintPenalties(torVec);
	}
	else if(pObjFun)
	{	
		for(unsigned k = 0; k < popSize; ++k)
		{
			setTOR(populationVec[k], *torTemplate);
			objFunValVec[k] = evaluateMultiObjective(torTemplate.get(), efF);
			addConstraintPenalties(torTemplate.get(), k);
			std::cout << "Done chromo #" << k << ", ofv: ";
			for(std::size_t kg = 0; kg < objFunValVec[k].first.size(); ++kg)
				std::cout << objFunValVec[k].first[kg] << " ";
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "Error: No objective function defined in GA optimizer" << std::endl;
		for(unsigned k = 0; k < popSize; ++k)
			objFunValVec[k].second = false;
	}
	if(numGoals > 1)
		fixOFVVecsForPareto();
}

void TopOptRealGA::addConstraintPenalties(const TopOptRep* const pTOR, unsigned k)
{
	std::pair<std::vector<double>, bool> cRes = evaluateMultiObjective(pTOR, efC);
	double addTerm = 0.;
	for(std::size_t kc = 0; kc < cRes.first.size(); ++kc)
		addTerm += pow(fabs(cRes.first[kc]), penaltyPower)*constraintPenalty;
	objFunValVec[k].first[0] += addTerm;
}

void TopOptRealGA::addConstraintPenalties(const std::vector<TopOptRep*>& torVec)
{
	assert(torVec.size() == popSize);
	for(unsigned k = 0; k < popSize; ++k)
		addConstraintPenalties(torVec[k], k);
}

void TopOptRealGA::doSelection()
{
	if(!noneValid())
	{
		if(numGoals == 1)
			GeneticOperators::Selection::tournamentSelection(populationVec, objFunValVec, ntourn);
		else if(numGoals == 2)
			GeneticOperators::Selection::paretoSelection2goals(populationVec, objFunValVec, goalWeights, sharingRadius);
	}
}

void debugPrint2D(const std::string& varName, std::list<double>& list1, const std::vector<std::size_t>& sizes)
{
	std::cout << varName << " = [";
	for(std::size_t k1 = 0; k1 < sizes[0]; ++k1)
	{
		for(std::size_t k2 = 0; k2 < sizes[1]; ++k2)
		{
			std::size_t curIndex = k2 + k1*sizes[1];
			std::list<double>::iterator lit1 = list1.begin();
			std::advance(lit1, curIndex);
			std::cout << *lit1 << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "];" << std::endl;
}

void TopOptRealGA::doCrossover()
{
	for(std::size_t kpop = 0; kpop < populationVec.size(); kpop+=2)
	{
		if(HelperNS::RandomGen::instance().coinFlip(crossRate) && (kpop + 1) < populationVec.size())
		{
			if(chromoSizes[1] == 1 || chromoSizes.size() == 3)
				GeneticOperators::Crossover::hybridize1D(populationVec[kpop], populationVec[kpop + 1]);
			else
				GeneticOperators::Crossover::hybridize2D(populationVec[kpop], populationVec[kpop + 1], chromoSizes);
			objFunValVec[kpop].second = false;
			objFunValVec[kpop + 1].second = false;
		}
	}
}

void TopOptRealGA::doMutation()
{
	for(std::size_t kpop = 0; kpop < populationVec.size(); ++kpop)
	{
		for(std::size_t kelem = 0; kelem < populationVec[kpop].size(); ++kelem)
		{
			if(HelperNS::RandomGen::instance().coinFlip(mutationRate))
			{
				if(chromoSizes[1] == 1 || chromoSizes.size() == 3)
					GeneticOperators::Mutation::standardMutation(populationVec[kpop], kelem, mutationRange);
				else
				{
					unsigned locMutRad = HelperNS::RandomGen::instance().randIntInRange((unsigned)0, curMutRad);
//					unsigned locMutRad = curMutRad;
					GeneticOperators::Mutation::nonlocalMutation2D(populationVec[kpop], chromoSizes, kelem, locMutRad, mutationRange);
				}
				objFunValVec[kpop].second = false;
			}
		}
	}
}

void TopOptRealGA::doElitism()
{
	if(!bestChromoVec.empty())
	{
		for(unsigned k = 0; k < numElite; ++k)
		{
			std::size_t krand = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, populationVec.size() - 1);
			std::size_t krand2 = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, bestChromoVec.size() - 1);
			populationVec[krand] = bestChromoVec[krand2];
			objFunValVec[krand].second = false;
			if(!objFunValVec[krand].first.empty())
			{
				objFunValVec[krand].first.resize(bestChromoOFVs[krand2].size());
				for(std::size_t k2 = 0; k2 < bestChromoOFVs[krand2].size(); ++k2)
					objFunValVec[krand].first[k2] = bestChromoOFVs[krand2][k2];
				objFunValVec[krand].second = true;
			}
		}
	}
}

void TopOptRealGA::boundsCheckChromos()
{
	for(std::size_t k = 0; k < populationVec.size(); ++k)
	{
		std::vector<double> realRep(populationVec[k].size());
		std::copy(populationVec[k].begin(), populationVec[k].end(), realRep.begin());
		torTemplate->boundsCheck(realRep);
		std::copy(realRep.begin(), realRep.end(), populationVec[k].begin());
	}
}

void TopOptRealGA::setTOR(const std::list<double>& x, TopOptRep& result) const
{
	std::vector<double> realRep(x.size());
	std::copy(x.begin(), x.end(), realRep.begin());
	result.setRealRep(realRep);
}

void TopOptRealGA::makeDiscreteChromos()
{
	for(std::size_t k = 0; k < populationVec.size(); ++k)
	{
		std::list<double>& curChromo = populationVec[k];
		HelperNS::greaterThanX locGTX(0.5);
		std::replace_if(curChromo.begin(), curChromo.end(), locGTX, 1.);
		HelperNS::lessThanX locLTX(0.51);
		std::replace_if(curChromo.begin(), curChromo.end(), locLTX, 0.);
	}
}

void TopOptRealGA::outputResults() const
{
	// Send best chromo to OutputHandlers
	if(!bestChromoVec.empty())
	{
		setTOR(bestChromoVec[0], *torTemplate);
		handleOutput(torTemplate.get());
	}
	// Output Pareto front for multi-objective problems
	if(numGoals > 1)
	{
		printBestChromos();
		printPopGoalSpace();
	}
}

void TopOptRealGA::printBestChromos() const
{
	std::ofstream outFile("bestChromos.m");
	for(std::size_t k = 0; k < bestChromoVec.size(); ++k)
	{
		setTOR(bestChromoVec[k], *torTemplate);
		std::unique_ptr<TOMesh> tmpMesh = torTemplate->getOutputMesh();
		OutputWriter::plotMeshMatlab(tmpMesh.get(), outFile);
		outFile << "title(num2str([";
		for(std::size_t kg = 0; kg < bestChromoOFVs[k].size(); ++kg)
			outFile << bestChromoOFVs[k][kg] << " ";
		outFile <<"]));" << std::endl;
	}
}

void TopOptRealGA::printPopGoalSpace() const
{
	std::ofstream outFile("paretoPop.m");
	outFile << "pop = [";
	for(std::size_t k = 0; k < objFunValVec.size(); ++k)
	{
		if(objFunValVec[k].second)
		{
			for(std::size_t k2 = 0; k2 < objFunValVec[k].first.size(); ++k2)
				outFile << objFunValVec[k].first[k2] << " ";
			outFile << std::endl;
		}
	}
	outFile << "];" << std::endl;
	outFile << "figure;" << std::endl;
	outFile << "plot(pop(:,1),pop(:,2),'b.');" << std::endl;
	outFile << "hold on" << std::endl;
	outFile << "best = [";
	for(std::size_t k = 0; k < bestChromoOFVs.size(); ++k)
	{
		for(std::size_t k2 = 0; k2 < bestChromoOFVs[k].size(); ++k2)
			outFile << bestChromoOFVs[k][k2] << " ";
		outFile << std::endl;
	}
	outFile << "];" << std::endl;
	outFile << "plot(best(:, 1), best(:, 2), 'ro');" << std::endl;
}

void TopOptRealGA::evaluateStatistics()
{
	if(!noneValid())
	{
		unsigned kstart = 0;
		while(!(objFunValVec[kstart].second && !objFunValVec[kstart].first.empty()) && kstart < popSize)
			kstart++;
		if(kstart >= popSize)
			return;
		numValid = 1;
		Real curBest, curWorst, sum = 0;
		unsigned bestChromoID = kstart, worstChromoID = kstart;
		curBest = objFunValVec[kstart].first[0];
		curWorst = curBest;
		sum = curBest;
		// find best & worst
		for(unsigned k = kstart + 1; k < popSize; k++)
		{
			if(objFunValVec[k].second && !objFunValVec[k].first.empty())
			{
				Real curVal = objFunValVec[k].first[0];
				if(curVal < curBest)
				{
					curBest = curVal;
					bestChromoID = k;
				}
				else if(curVal > curWorst)
				{
					curWorst = curVal;
					worstChromoID = k;
				}
				sum += curVal;
				numValid++;
			}
		}
		popBest = curBest;
		if(popBest < bestEver || bestChromoVec.empty())
		{
			bestEver = popBest;
			bestChromoVec.resize(1);
			bestChromoVec[0] = populationVec[bestChromoID];
			bestChromoOFVs.resize(1);
			bestChromoOFVs[0] = objFunValVec[bestChromoID].first;
		}
		popAvg = sum/(double)numValid;
		// Compute standard deviation
		Real stdDev = 0;
		for(unsigned k = kstart; k < popSize; k++)
		{
			if(objFunValVec[k].second && !objFunValVec[k].first.empty())
			{
				Real curVal = objFunValVec[k].first[0];
				stdDev += (curVal - popAvg)*(curVal - popAvg);
			}
		}
		popStdDev = sqrt(stdDev/(Real)numValid);

		if(numGoals > 1)
			updateBestChromosPareto();
	}
	else
		numValid = 0;
}

void TopOptRealGA::updateBestChromosPareto()
{
	// Append old best chromos with current pop
	std::vector<std::pair<std::vector<double>, bool> > objFunOfPopAndBestVec;
	combineOFVVecs(objFunOfPopAndBestVec);
	// Compute rank of combined vectors
	std::vector<std::size_t> rank;
	GeneticOperators::Selection::computeParetoRank(rank, objFunOfPopAndBestVec);
	std::vector<std::list<double> > prevBestChromoVec = bestChromoVec;
	std::vector<std::vector<double> > prevBestOFVs = bestChromoOFVs;
	bestChromoVec.clear();
	bestChromoOFVs.clear();
	for(std::size_t k = 0; k < populationVec.size(); ++k)
	{
		if(rank[k] == 1)
		{
			bestChromoVec.push_back(populationVec[k]);
			bestChromoOFVs.push_back(objFunValVec[k].first);
		}
	}
	std::size_t start = populationVec.size();
	for(std::size_t k = start; k < rank.size(); ++k)
	{
		if(rank[k] == 1)
		{
			bestChromoVec.push_back(prevBestChromoVec[k - start]);
			bestChromoOFVs.push_back(prevBestOFVs[k - start]);
		}
	}
	removeDuplicateBestChromos();
}

void TopOptRealGA::combineOFVVecs(std::vector<std::pair<std::vector<double>, bool> >& outVec) const
{
	// FIrst generate a bestOFVvec compatible with objFunValVec
	outVec.clear();
	outVec.reserve(objFunValVec.size() + bestChromoOFVs.size());
	outVec = objFunValVec;
	for(std::size_t k = 0; k < bestChromoOFVs.size(); ++k)
	{
		std::pair<std::vector<double>, bool> curElem(bestChromoOFVs[k], true);
		outVec.push_back(curElem);
	}
}

void TopOptRealGA::removeDuplicateBestChromos()
{
	for(std::size_t k = 0; k < bestChromoOFVs.size(); ++k)
	{
		for(std::size_t k2 = k + 1; k2 < bestChromoOFVs.size(); ++k2)
		{
			if(bestChromoOFVs[k] == bestChromoOFVs[k2])
			{
				bestChromoOFVs.erase(bestChromoOFVs.begin() + k2);
				bestChromoVec.erase(bestChromoVec.begin() + k2);
			}
		}
	}
}

void TopOptRealGA::printStatistics(std::ofstream& outFile) const
{
	if(numValid > 0)
	{
		std::cout << "Number of valid chromos: " << numValid << std::endl;
		std::cout << "Population average: " << popAvg << std::endl;
		std::cout << "Population standard deviation: " << popStdDev << std::endl;
		std::cout << "Generation best: " << popBest << std::endl;
		std::cout << "Best ever: " << bestEver << std::endl;
	}
	else
		std::cout << "No valid chromosomes were evaluated." << std::endl;
	// file output
	if(numValid > 0)
	{
		outFile << numValid << " ";
		outFile << bestEver << " ";
		outFile << popBest << " ";
		outFile << popAvg << " ";
		outFile << popStdDev << std::endl;
	}
	else
		outFile << "No valid chromosomes!" << std::endl;
}

bool TopOptRealGA::noneValid() const
{
	for(std::size_t k = 0; k < objFunValVec.size(); ++k)
	{
		if(objFunValVec[k].second)
			return false;
	}
	return true;
}

void TopOptRealGA::fixOFVVecsForPareto()
{
	for(std::size_t k = 0; k < objFunValVec.size(); ++k)
	{
		if(objFunValVec[k].first.size() < numGoals)
			objFunValVec[k].second = false;
	}
}
}

