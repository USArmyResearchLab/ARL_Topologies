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

#ifndef TORFACTORY_H
#define TORFACTORY_H

#include <memory>
#include "inputloader.h"

namespace Topologies{
class TopOptRep;
class TopOptObjFun;
class MPIHandler;

//! Namespace for generating TopOptRep derived classes from input files and data vectors
namespace TopOptRepFactory
{
	//! Returns a TopOptRep object as specified by a RepNodeInfo object
	std::unique_ptr<TopOptRep> createTopOptRep(const InputLoader::RepNodeInfo& rni);
	//! Returns a TopOptRep object as specified by a TORType and an OptNodeInfo object
  std::unique_ptr<TopOptRep> createTopOptRep(TORType inTORT, const InputLoader::OptNodeInfo& oni);
	//! Returns a TopOptRep object as specified by the class name in `torName` and the given data
	std::unique_ptr<TopOptRep> createTopOptRep(TORType inTORT, const std::vector<std::vector<double> >& realRep,
                                             const std::vector<std::vector<int> >& discreteRep);
}
}
#endif

