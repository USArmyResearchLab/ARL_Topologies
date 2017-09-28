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

#ifndef TOPOPTFACTORY_H
#define TOPOPTFACTORY_H

#include <memory>
#include "inputloaderopt.h"

namespace Topologies{
class TopOpt;
class TopOptObjFun;
class MPIHandler;
class OutputHandler;

//! Namespace containing factory functions for creating TopOpt objects
namespace TopOptFactory
{
	//! Returns a TopOpt object with specified input parameters as given in the file inputFileName
	std::unique_ptr<TopOpt> createTopOpt(const InputLoader::OptNodeInfo& oni, TopOptObjFun* pObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inpMPIH = nullptr);
}
}
#endif

