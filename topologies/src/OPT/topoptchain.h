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

#ifndef TOPOPTCHAIN_H
#define TOPOPTCHAIN_H

#include "topopt.h"
#include <memory>
#include <vector>
#include <string>

namespace Topologies{
class TopOptRep;
class TopOptObjFun;
class MPIHandler;
class OutputHandler;

//! An optimization class for strining together a set of optimizers and other operations.
/*! This class is used for implementing continuation methods and refinement.  A set of optimizers
 *  each with different input parameters can be set up and run.  This is mostly useful for 
 *  performing continuation on TopOptRep values such as penalty and Heaviside beta.  Refinement
 *  can also be performed for PixelRep and VoxelRep.
 */
class TopOptChain : public TopOpt
{
public:
	TopOptChain(const InputLoader::TOOChain& inputData, TopOptObjFun* inpObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inpMPIH = nullptr);
	virtual ~TopOptChain();
	TopOptChain(const TopOptChain& copy);
  TopOptChain(TopOptChain&& copy) {swap(copy);}
  TopOptChain& operator=(TopOptChain rhs) {swap(rhs); return *this;}
  void swap(TopOptChain& arg);

	virtual std::unique_ptr<TopOptRep> optimize(const TopOptRep& initialGuess);

private:
	void setupTOO(const InputLoader::OptNodeInfo& oni);
	void refineTOR(TopOptRep* torToRefine);
	std::unique_ptr<TopOptRep> convertTOR(const TopOptRep& torToConvert, const InputLoader::OptNodeInfo& oni);
	std::unique_ptr<TopOptRep> continuationTOR(const TopOptRep& torToConvert, const InputLoader::OptNodeInfo& oni);

	std::vector<InputLoader::OptNodeInfo> oniVec;
	std::unique_ptr<TopOpt> toOptimizer;
};

inline
TopOptChain::TopOptChain(const TopOptChain& copy) :
	oniVec(copy.oniVec),
	toOptimizer(nullptr) // Note: this shouldn't matter as calling optimize will set this ptr
{}

inline
void TopOptChain::swap(TopOptChain& arg)
{
	oniVec.swap(arg.oniVec);
	toOptimizer.swap(arg.toOptimizer);
}
}
#endif

