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

#include "geneticoperators.h"
#include "helper.h"

namespace Topologies{
namespace GeneticOperators
{
	namespace Crossover
	{
		void hybridize1D(std::list<double>& chromo1, std::list<double>& chromo2)
		{
		  if(chromo1.size() != chromo2.size() || chromo1.empty() || chromo2.empty())
			{
				std::cout << "Error in input to hybridize1D!" << std::endl;
				return;
			}
			// Choose a random point in the lists
			std::size_t krand = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, chromo1.size() - 1);
			std::list<double>::iterator it1 = chromo1.begin(), it2 = chromo2.begin();
			std::advance(it1, krand);
			std::advance(it2, krand);
			// Hybridize the numbers at these points
			double d1 = *it1, d2 = *it2;
			double r = HelperNS::RandomGen::instance().randRealInRange<double>();
			*it1 = r*d1 + (1. - r)*d2;
			*it2 = r*d2 + (1. - r)*d1;
			// Splice lists
			chromo1.splice(it1, chromo2, it2, chromo2.end());
			chromo2.splice(chromo2.end(), chromo1, it1, chromo1.end());
		}

		void hybridize2D(std::list<double>& chromo1, std::list<double>& chromo2, const std::vector<std::size_t>& sizes)
		{
			if(chromo1.size() != chromo2.size() || chromo1.empty() || chromo2.empty() || sizes.empty())
			{
				std::cout << "Error in input to hybridize2D!" << std::endl;
				return;
			}
			// Pick rows or columns
			if(HelperNS::RandomGen::instance().coinFlip())
  		{
				// Hybridize column
				double r = HelperNS::RandomGen::instance().randRealInRange<double>();
				std::size_t krand = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, sizes[1] - 1);
				for(std::size_t k1 = 0; k1 < sizes[0]; ++k1)
				{
					std::size_t curIndex = krand + k1*sizes[1];
					std::list<double>::iterator lit1 = chromo1.begin();
					std::advance(lit1, curIndex);
					std::list<double>::iterator lit2 = chromo2.begin();
					std::advance(lit2, curIndex);
					double d1 = *lit1, d2 = *lit2;
					*lit1 = r*d1 + (1. - r)*d2;
					*lit2 = r*d2 + (1. - r)*d1;
    		}
				// Swap columns
				for(std::size_t k1 = 0; k1 < sizes[0]; ++k1)
				{
					for(std::size_t k2 = krand; k2 < sizes[1]; ++k2)
					{
						std::size_t curIndex = k2 + k1*sizes[1];
						std::list<double>::iterator lit1 = chromo1.begin();
						std::advance(lit1, curIndex);
						std::list<double>::iterator lit2 = chromo2.begin();
						std::advance(lit2, curIndex);
						std::swap(*lit1, *lit2);
					}
				}
			}
			else
			{
				// Hybridize row
				double r = HelperNS::RandomGen::instance().randRealInRange<double>();
				std::size_t krand = HelperNS::RandomGen::instance().randIntInRange((std::size_t)0, sizes[0] - 1);
				for(std::size_t k2 = 0; k2 < sizes[1]; ++k2)
				{
					std::size_t curIndex = k2 + krand*sizes[1];
					std::list<double>::iterator lit1 = chromo1.begin();
					std::advance(lit1, curIndex);
					std::list<double>::iterator lit2 = chromo2.begin();
					std::advance(lit2, curIndex);
					double d1 = *lit1, d2 = *lit2;
					*lit1 = r*d1 + (1. - r)*d2;
					*lit2 = r*d2 + (1. - r)*d1;
				}
				// Swap rows
				for(std::size_t k1 = krand; k1 < sizes[0]; ++k1)
				{
					for(std::size_t k2 = 0; k2 < sizes[1]; ++k2)
					{
						std::size_t curIndex = k2 + k1*sizes[1];
						std::list<double>::iterator lit1 = chromo1.begin();
						std::advance(lit1, curIndex);
						std::list<double>::iterator lit2 = chromo2.begin();
						std::advance(lit2, curIndex);
						std::swap(*lit1, *lit2);
					}
				}
			}
		}

		void hybridize3D(std::list<double>& chromo1, std::list<double>& chromo2, const std::vector<std::size_t>& sizes)
		{
			std::cout << "Not yet implemented" << std::endl;
			abort();
		}
	}

	namespace Mutation
	{
		void standardMutation(std::list<double>& chromo, std::size_t kelem, double range)
		{
			double rval = HelperNS::RandomGen::instance().randRealInRange(-range, range);
			std::list<double>::iterator kit = chromo.begin();
			std::advance(kit, kelem);
			double curval = *kit;
			*kit = curval + rval;
		}

    void nonlocalMutation2D(std::list<double>& chromo, const std::vector<std::size_t>& sizes, std::size_t kelem, 
														unsigned radius, double range)
		{
			unsigned rad2 = radius*radius;
			std::size_t e1 = kelem/sizes[1];
			std::size_t e2 = kelem%sizes[1];
			double rval = HelperNS::RandomGen::instance().randRealInRange(-range, range);
			std::list<double>::iterator kit = chromo.begin();
			if(chromo.size() != sizes[0]*sizes[1])
			{
				std::cout << "Error, chromo size doesn't match TOR template size in GA mutation" << std::endl;
				return;
			}
			for(std::size_t k1 = 0; k1 < sizes[0]; ++k1)
			{
				for(std::size_t k2 = 0; k2 < sizes[1]; ++k2)
				{
					std::size_t d1 = HelperNS::absdiff(k1, e1);
					std::size_t d2 = HelperNS::absdiff(k2, e2);
					if((d1*d1 + d2*d2) <= rad2)
					{
						double curval = *kit;
						*kit = curval + rval;
					}
      		++kit;
				}
			}
		}
	}
}
}
