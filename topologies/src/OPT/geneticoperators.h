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

#ifndef GENETICOPERATORS_H
#define GENETICOPERATORS_H

#include <vector>
#include <list>

namespace Topologies{
//! Namespace that contains functions for TopOptGA
namespace GeneticOperators
{
	//! Namespace that contains crossover opterators for TopOptGA
	namespace Crossover
	{
		//! 1D, single point crossover for chromo1 and chromo2
		/*! This function treats the input lists as 1D and randomly chooses a crossover point.
		 *  The values at the crossover point are hybridized and the values past that point are swapped.
		 */
		void hybridize1D(std::list<double>& chromo1, std::list<double>& chromo2);
		//! 2D, single point crossover for chromo1 and chromo2
		/*! This function treats the input lists as 2D with sizes given in sizes.  A random row or column
		 *  is chosen, the values along it are hybridized and the values past that row/column are swapped.
		 */
		void hybridize2D(std::list<double>& chromo1, std::list<double>& chromo2, const std::vector<std::size_t>& sizes);
		//! 3D, single point crossover for chromo1 and chromo2
		/*! This function treats the input lists as 2D with sizes given in sizes.  A random plane is chosen,
		 *  the values along it are hybridized and the values past that row/column are swapped.
		 */
		void hybridize3D(std::list<double>& chromo1, std::list<double>& chromo2, const std::vector<std::size_t>& sizes);
		//! Uniform crossover implementation
		/*! This implementation swaps a random subset of values in chromo1 and chromo2 and blends a random subset of the
		 *  swapped values.  The size of the subsets is chosen randomly.
		 */
		template<typename container>
		void hybridizeUniform(container& chromo1, container& chromo2);
	}

	//! Namespace that contains mutation opterators for TopOptGA
	namespace Mutation
	{
		//! Mutation function that adds a random value with maximum value range to the kelem element of chromo
		template<typename container>
		void standardMutation(container& chromo, std::size_t kelem, double range);
		//! Similar to standardMutation, but this function also mutates the surrounding values (assuming a 2d chromo) within radius distance
		template<typename container>
		void nonlocalMutation2D(container& chromo, const std::vector<std::size_t>& sizes, std::size_t kelem, 
														unsigned radius, double range);
		//! Similar to standardMutation, but this function also mutates the surrounding values (assuming a 3d chromo) within radius distance
		template<typename container>
		void nonlocalMutation3D(container& chromo, const std::vector<std::size_t>& sizes, std::size_t kelem, 
														unsigned radius, double range);
	}
}
}
#endif

