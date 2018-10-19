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

#ifndef INPUTLOADEROPT_H
#define INPUTLOADEROPT_H

#include "inputloader.h"

namespace Topologies{
namespace InputLoader
{
	//! Class to handle loading of TopOptChain
	class TOOChain : public InputParser
	{
	public:
		TOOChain() {}
		virtual ~TOOChain() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns an iterator to the first OptNodeInfo object
		std::vector<OptNodeInfo>::const_iterator oniBegin() const {return oniVec.begin();}
		//! Returns the past-the-end iterator of the OptNodeInfo object vector
		std::vector<OptNodeInfo>::const_iterator oniEnd() const {return oniVec.end();}
	private:
		std::vector<OptNodeInfo> oniVec;
	};

	//! Class to handle loading of TopOptOC, TopOptNLOpt, and TopOptNLOptUnc
	class TOOGeneric : public InputParser
	{
	public:
		TOOGeneric(const std::string& inOptName) : optimizerName(inOptName) {}
		virtual ~TOOGeneric() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns the name of the optimizer
		const std::string& getOptimizerName() const {return optimizerName;}
		//! Returns the gradient filter size
		double getFilterSize() const {return filterSize;}
		//! Returns the constraint penalty
		double getConstraintPenalty() const {return constraintPenalty;}
		//! Returns the power to use on the constraint penalty
		double getPenaltyPower() const {return penaltyPower;}
		//! Returns the step size
		double getStepSize() const {return stepSize;}
		//! Returns the maximum step size
		double getMaxStep() const {return maxStep;}
		//! Returns the stop tolerance
		double getStopTol() const {return stopTol;}
		//! Returns the maximum number of iterations
		unsigned getMaxIters() const {return maxIters;}
		//! Returns the number of constraints
		unsigned getNumConstraints() const {return numConstraints;}
	private:
		const std::string optimizerName;
		double filterSize = 0., constraintPenalty = 0., penaltyPower = 1.;
		double stepSize = 1., maxStep = 1.;
		double stopTol = 0.01;
		unsigned maxIters = 100, numConstraints = 1;
	};

	//! Class to handle loading of TopOptGA and TopOptPGA
	class TOOGA : public InputParser
	{
	public:
		TOOGA(const std::string& inOptName) : isPareto(inOptName == "pga") {}
		virtual ~TOOGA() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns whether or not this is a Pareto optimizer
		bool getIsPareto() const {return isPareto;}
		//! Returns the constraint penalty
		double getConstraintPenalty() const {return constraintPenalty;}
		//! Returns the power applied to the constraint penalty
		double getPenaltyPower() const {return penaltyPower;}
		//! Returns the mutation range
		double getMutationRange() const {return mutationRange;}
		//! Returns the crossover rate
		double getCrossoverRate() const {return crossRate;}
		//! Returns the mutation rate
		double getMutationRate() const {return mutationRate;}
		//! Returns the population size
		unsigned getPopSize() const {return popSize;}
		//! Returns the number of generations
		unsigned getNumGens() const {return numGens;}
		//! Returns the number of chromosomes in a tournament for selection
		unsigned getNumTourn() const {return ntourn;}
		//! Returns the number of chromosomes to replace for elitism
		unsigned getNumElite() const {return numElite;}
		//! Returns the mutation radius
		unsigned getMutationRadius() const {return  mutationRadius;}
		//! Returns the number of goals for optimization
		unsigned getNumGoals() const {return numGoals;}
		//! Returns the sharing radius (for Pareto)
		double getSharingRadius() const {return sharingRadius;}
		//! Returns the weighting values for each goal (for Pareto)
		const std::vector<double>&  getGoalWeights() const {return goalWeights;}
	private:
		const bool isPareto;
		double constraintPenalty = 0., penaltyPower = 1., mutationRange = 1.;
		double crossRate = 1., mutationRate = 0.;
		unsigned popSize = 2, numGens = 1, ntourn = 2, numElite = 1, mutationRadius = 1;
		unsigned numGoals = 1;
		double sharingRadius = 0.1;
		std::vector<double> goalWeights;
	};
}
}
#endif

