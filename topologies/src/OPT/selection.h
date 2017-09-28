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

#ifndef SELECTION_H
#define SELECTION_H

#include "geneticoperators.h"
#include <iostream>

namespace Topologies{
namespace GeneticOperators
{
	//! A collection of functions for GA selection, used in TopOptGA
	namespace Selection
	{
		//! Tournament selection chooses ntourn chromosomes, and the best chromosome of the set is copied to the new popoulation
		template <typename T>
		void tournamentSelection(std::vector<T>& populationVec, std::vector<std::pair<std::vector<double>, bool> >& objFunValVec, 
														unsigned ntourn)
		{
			std::vector<std::pair<std::vector<double>, bool> > copyOFVVec = objFunValVec;
			std::vector<T> copyPopVec = populationVec;
			std::size_t npop = populationVec.size();
			for(std::size_t kpop = 0; kpop < npop; ++kpop)
			{
				std::size_t krand = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, npop - 1);
				while(!copyOFVVec[krand].second) // find valid chromosome
					krand = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, npop - 1);
				std::size_t curBest = krand;
				unsigned numChecked = 0;
				while(numChecked < ntourn)
				{
					krand = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, npop - 1);
					if(copyOFVVec[krand].second)
					{
						numChecked++;
						if(copyOFVVec[krand].first < copyOFVVec[curBest].first)
							curBest = krand;
					}
				}
				// Copy chromosome to new population
				populationVec[kpop] = copyPopVec[curBest];
				objFunValVec[kpop] = copyOFVVec[curBest];
			}
		}

		namespace
		{
			typedef std::pair<std::vector<double>, bool> OFV;
			// Anonymous helper functions for paretoSelection2goals
			bool allChromosRanked(const std::vector<std::size_t>& rank,	const std::vector<OFV>& objFunValVec)
			{
				for(int kpop = 0; kpop < rank.size(); kpop++)
				{
					if(rank[kpop] == 0 && objFunValVec[kpop].second)
						return false;
				}
				return true;
			}

			void paretoRankPopulation(std::vector<std::size_t>& rank, const std::vector<OFV>& objFunValVec, std::size_t curRank)
			{
				std::size_t popSize = objFunValVec.size();
				for(int kpop1 = 0; kpop1 < popSize; kpop1++)
				{
					const OFV& ofv1 = objFunValVec[kpop1];
					bool dominated = !ofv1.second;
					if(rank[kpop1] == 0 && ofv1.second)
					{
						for(std::size_t kpop2 = 0; kpop2 < popSize && !dominated; kpop2++)
						{
							if(rank[kpop2] == 0 || rank[kpop2] == curRank)
							{
								const OFV& ofv2 = objFunValVec[kpop2];
								if((ofv2.first[0] < ofv1.first[0] && ofv2.first[1] < ofv1.first[1]) ||
								   (ofv2.first[0] <= ofv1.first[0] && ofv2.first[1] < ofv1.first[1]) ||
								   (ofv2.first[0] < ofv1.first[0] && ofv2.first[1] <= ofv1.first[1]))
								{
									dominated = true;
								}
							}
						}
						if(!dominated)
							rank[kpop1] = curRank;
					}
				}
			}

			void computeSharingVals(std::vector<double>& sharingVals, const std::vector<OFV>& objFunValVec, 
				const std::vector<double>& paretoGoalWeights, double sharingRadius)
			{
				Real sum = 0;
				std::size_t popSize = objFunValVec.size();
				sharingVals.resize(popSize);
				for(int kpop1 = 0; kpop1 < popSize; kpop1++)
				{
					sharingVals[kpop1] = 1.;
					const OFV& ofv1 = objFunValVec[kpop1];
					double f11 = ofv1.first[0]/paretoGoalWeights[0], f12 = ofv1.first[1]/paretoGoalWeights[1];
					for(int kpop2 = 0; kpop2 < popSize; kpop2++)
					{
						const OFV& ofv2 = objFunValVec[kpop2];
						Real f21 = ofv2.first[0]/paretoGoalWeights[0], f22 = ofv2.first[1]/paretoGoalWeights[1];
						Real d = sqrt(((f11 - f21)*(f11 - f21)) + (f12 - f22)*(f12 - f22));
						if(d < sharingRadius && ofv2.second)
							sharingVals[kpop1] += 1. - d/sharingRadius;
					}
					sum += sharingVals[kpop1];
				}
			}

			void computeParetoOFV(const std::vector<Real>& sharingVals, const std::vector<std::size_t>& rank, 
				std::vector<double>& ofvVec, std::size_t maxRank, bool useSharing)
			{
				Real curF = 1., minOFV = -1;
				std::size_t popSize = sharingVals.size();
				ofvVec.resize(popSize);
				for(int krank = 1; krank < maxRank; krank++)
				{
					for(int kpop = 0; kpop < popSize; kpop++)
					{
						if(rank[kpop] == krank)
						{
							if(!useSharing)
								ofvVec[kpop] = 1./static_cast<Real>(krank);
							else
								ofvVec[kpop] = curF/sharingVals[kpop];
							if(minOFV < 0 || (ofvVec[kpop] < minOFV))
								minOFV = ofvVec[kpop];
						}
					}
					curF = minOFV;
					minOFV = -1;
				}
			}

			template <typename T>
			void rouletteSelection(const std::vector<double> paretoOFVVec, std::vector<T>& populationVec)
			{
				std::vector<T> copyPop = populationVec;
				Real ofvSum = 0;
				std::size_t popSize = paretoOFVVec.size();
				for(int kpop = 0; kpop < popSize; kpop++)
					ofvSum += paretoOFVVec[kpop];
				for(std::size_t kpop = 0; kpop < popSize; kpop++)
				{
					Real tmpURV = HelperNS::RandomGen::instance().randRealInRange((Real)0., ofvSum);
					Real tmpOFV = paretoOFVVec[0];
					int kselect = 0;
					while(tmpOFV < tmpURV && kselect < popSize)
					{
						kselect++;
						tmpOFV += paretoOFVVec[kselect];
					}
					populationVec[kpop] = copyPop[kselect];
				}
			}
		}

		//! Compute the Pareto rank of each chromosome with specified objective function values
		/*! The Pareto rank is a value used for giving a sense of fitness of a chromosome.
		 *  It essentially gives the number of sets of chromosomes that dominate it (plus 1).  In
		 *  other words, Pareto optimal chromosomes are rank 1 as they are non-dominated.  These
		 *  Pareto optimal chromos are removed and the Pareto optimal chromos from the remaining set
		 *  are determined, which are the rank 2 chromos.  The process is repeated until all chromos have
		 *  a rank.
		 */
		std::size_t computeParetoRank(std::vector<std::size_t>& rank, 
			const std::vector<std::pair<std::vector<double>, bool> >& objFunValVec)
		{
			std::size_t npop = objFunValVec.size();
			rank.resize(npop);
			std::size_t curRank = 1;
			for(std::size_t kpop = 0; kpop < npop; kpop++)
				rank[kpop] = 0;
			while(!allChromosRanked(rank, objFunValVec) && curRank < 1000)
			{
				paretoRankPopulation(rank, objFunValVec, curRank);
				curRank++;
			}
			// Check if any remain unranked
			for(std::size_t kpop = 0; kpop < npop; kpop++)
				if(rank[kpop] == 0)
					rank[kpop] = curRank;
			// set invalid chromos to last rank
			for(int kpop = 0; kpop < npop; kpop++)
				if(!objFunValVec[kpop].second)
					rank[kpop] = curRank;
			curRank++;
			return curRank;
		}

		//! A selection routine for Pareto optimization of 2 goals
		/*! This method first computes the Pareto rank of all chromos, next computes
		 *  sharing values, which give a measure of how crowded a chromosome is in goal-space, and then
		 *  computes a single objective value based on rank and crowding (lower rank and lower crowding give
		 *  higher values).  Roulette wheel selection is used on this value.
		 */
		template <typename T>
		void paretoSelection2goals(std::vector<T>& populationVec, 
			std::vector<std::pair<std::vector<double>, bool> >& objFunValVec, 
			const std::vector<double>& paretoGoalWeights, double sharingRadius)
		{
			// Check validity
			if(objFunValVec.empty())
				return;
			std::size_t npop = objFunValVec.size();
			if(objFunValVec[0].first.size() != 2)
			{
				std::cout << "Error: paretoSelection2goals only valid for 2 goals, found " << objFunValVec[0].first.size() << std::endl;
				return;
			}
			// Rank population
			std::vector<std::size_t> rank(npop);
			std::size_t maxRank = computeParetoRank(rank, objFunValVec);
			// compute sharing values
			std::vector<double> sharingVals(npop);
			computeSharingVals(sharingVals, objFunValVec, paretoGoalWeights, sharingRadius);
			// compute new objective function values
			std::vector<double> paretoOFVVec;
			computeParetoOFV(sharingVals, rank, paretoOFVVec, maxRank, true);
			// perform roulette wheel selection
			rouletteSelection(paretoOFVVec, populationVec);
		}
	}
}
}
#endif

