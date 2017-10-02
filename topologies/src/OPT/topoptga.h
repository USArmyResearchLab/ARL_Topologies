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

#ifndef TOPOPTGA_H
#define TOPOPTGA_H

#include "topopt.h"
#include "topoptrep.h"
#include <memory>
#include <vector>
#include <fstream>
#include <list>

namespace Topologies{
class TopOptObjFun;
class MPIHandler;
class OutputHandler;

//! A basic real-coded GA that also handles multi-objective optimization
class TopOptRealGA : public TopOpt
{
public:
	TopOptRealGA(const InputLoader::TOOGA& inputData, TopOptObjFun* inpObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inpMPIH = nullptr);
	virtual ~TopOptRealGA();
	TopOptRealGA(const TopOptRealGA& copy);
	TopOptRealGA(TopOptRealGA&& copy) {swap(copy);}
	TopOptRealGA& operator=(TopOptRealGA rhs) {swap(rhs); return *this;}
	void swap(TopOptRealGA& arg);

	virtual std::unique_ptr<TopOptRep> optimize(const TopOptRep& initialGuess);

private:
	void initializePopulation(TopOptRep& initialGuess);
	void evaluatePop();
	void addVolFracPenalties();
	void doSelection();
	void doElitism();
	void doCrossover();
	void doMutation();
	void setTOR(const std::list<double>& x, TopOptRep& result) const;
	void outputResults() const;
	void printBestChromos() const;
	bool noneValid() const;
	void printStatistics(std::ofstream& outFile) const;
	void evaluateStatistics();
	void fixOFVVecsForPareto();
	void printPopGoalSpace() const;
	void combineOFVVecs(std::vector<std::pair<std::vector<double>, bool> >& outVec) const;
	void removeDuplicateBestChromos();
	void updateBestChromosPareto();
	void boundsCheckChromos();
	void makeDiscreteChromos();
	void addConstraintPenalties(const std::vector<TopOptRep*>& torVec);
	void addConstraintPenalties(const TopOptRep* const pTOR, unsigned k);
	// Data
	double mutationRange;
	double constraintPenalty, penaltyPower;
	double crossRate, mutationRate;
	bool useDiscreteGA;
	unsigned popSize, numGens, ntourn, numElite;
	unsigned maxMutRad, curMutRad;
	mutable std::unique_ptr<TopOptRep> torTemplate;
	std::vector<std::size_t> chromoSizes;
	// Pareto stuff
	unsigned numGoals;
	std::vector<double> goalWeights;
	double sharingRadius;
	// population
	std::vector<std::list<double> > populationVec;
	std::vector<std::pair<std::vector<double>, bool> > objFunValVec;
	std::vector<std::list<double> > bestChromoVec;
	std::vector<std::vector<double> > bestChromoOFVs;
	// Stats
	double popAvg, popBest, bestEver, popStdDev;
	unsigned numValid;
};

inline
TopOptRealGA::TopOptRealGA(const TopOptRealGA& copy) :
	TopOpt(copy),
	mutationRange(copy.mutationRange),
	constraintPenalty(copy.constraintPenalty),
	penaltyPower(copy.penaltyPower),
	mutationRate(copy.mutationRate),
	crossRate(copy.crossRate),
	useDiscreteGA(copy.useDiscreteGA),
	popSize(copy.popSize),
	numGens(copy.numGens),
	ntourn(copy.ntourn),
	numElite(copy.numElite),
	maxMutRad(copy.maxMutRad),
	curMutRad(copy.curMutRad),
	torTemplate(copy.torTemplate->clone()),
	chromoSizes(copy.chromoSizes),
	numGoals(copy.numGoals),
	goalWeights(copy.goalWeights),
	sharingRadius(copy.sharingRadius),
	populationVec(copy.populationVec),
	objFunValVec(copy.objFunValVec),
	bestChromoVec(copy.bestChromoVec),
	bestChromoOFVs(copy.bestChromoOFVs),
	popAvg(copy.popAvg),
	popBest(copy.popBest),
	bestEver(copy.bestEver),
	popStdDev(copy.popStdDev),
	numValid(copy.numValid)
{
}

inline
void TopOptRealGA::swap(TopOptRealGA& arg)
{
	TopOpt::swap(arg);
	std::swap(mutationRange, arg.mutationRange);
	std::swap(constraintPenalty, arg.constraintPenalty);
	std::swap(penaltyPower, arg.penaltyPower);
	std::swap(crossRate, arg.crossRate);
	std::swap(mutationRate, arg.mutationRate);
	std::swap(useDiscreteGA, arg.useDiscreteGA);
	std::swap(popSize, arg.popSize);
	std::swap(numGens, arg.numGens);
	std::swap(ntourn, arg.ntourn);
	std::swap(numElite, arg.numElite);
	std::swap(maxMutRad, arg.maxMutRad);
	std::swap(curMutRad, arg.curMutRad);
	torTemplate.swap(arg.torTemplate);
	chromoSizes.swap(arg.chromoSizes);
	std::swap(numGoals, arg.numGoals);
  std::swap(goalWeights, arg.goalWeights);
  std::swap(sharingRadius, arg.sharingRadius);
  populationVec.swap(arg.populationVec);
  objFunValVec.swap(arg.objFunValVec);
  bestChromoVec.swap(arg.bestChromoVec);
  bestChromoOFVs.swap(arg.bestChromoOFVs);
  std::swap(popAvg, arg.popAvg);
	std::swap(popBest, arg.popBest);
	std::swap(bestEver, arg.bestEver);
	std::swap(popStdDev, arg.popStdDev);
	std::swap(numValid, arg.numValid);
}
}
#endif


