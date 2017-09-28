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

#include "topoptfactory.h"
#include "topopt.h"
#include "topoptoc.h"
#include "topoptnlopt.h"
#include "topoptga.h"
#include "topoptgd.h"
#include "topoptchain.h"
#include "inputloader.h"

namespace Topologies{
namespace TopOptFactory
{
	std::unique_ptr<TopOpt> createTopOpt(const InputLoader::OptNodeInfo& oni, TopOptObjFun* pObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inpMPIH)
	{
		std::unique_ptr<TopOpt> toOptimizer;
		if(oni.getType() == tootOC)
		{
			InputLoader::TOOGeneric tooParser(oni.getTypeName());
			tooParser.parseNode(oni);
			toOptimizer = std::unique_ptr<TopOpt>(new TopOptOC(tooParser, pObjFun, inOutVec, inpMPIH));
		}
		else if(oni.getType() == tootMMA)
		{
			InputLoader::TOOGeneric tooParser(oni.getTypeName());
			tooParser.parseNode(oni);
			toOptimizer = std::unique_ptr<TopOpt>(new TopOptNLOpt(oni.getType(), tooParser, pObjFun, inOutVec, inpMPIH));
		}
		else if(oni.getType() == tootBFGS)
		{
			InputLoader::TOOGeneric tooParser(oni.getTypeName());
			tooParser.parseNode(oni);
			toOptimizer = std::unique_ptr<TopOpt>(new TopOptNLOpt(oni.getType(), tooParser, pObjFun, inOutVec, inpMPIH));
		}
		else if(oni.getType() == tootGA || oni.getType() == tootPGA)
		{
			InputLoader::TOOGA tooParser(oni.getTypeName());
			tooParser.parseNode(oni);
			toOptimizer = std::unique_ptr<TopOpt>(new TopOptRealGA(tooParser, pObjFun, inOutVec, inpMPIH));
		}
		else if(oni.getType() == tootGD)
		{
			InputLoader::TOOGeneric tooParser(oni.getTypeName());
			tooParser.parseNode(oni);
			toOptimizer = std::unique_ptr<TopOpt>(new TopOptGD(tooParser, pObjFun, inOutVec, inpMPIH));
		}
		else if(oni.getType() == tootChain)
		{
			InputLoader::TOOChain tooParser;
			tooParser.parseNode(oni);
			toOptimizer = std::unique_ptr<TopOpt>(new TopOptChain(tooParser, pObjFun, inOutVec, inpMPIH));
		}
		else if(oni.getType() == tootRefine)
		{
			std::cout << "Received optimizer REFINE, which should only be used with CHAIN, aborting" << std::endl;
			abort();
		}
		else if(oni.getType() == tootConvert)
		{
			std::cout << "Received optimizer CONVERT, which should only be used with CHAIN, aborting" << std::endl;
			abort();
		}
		else
		{
			std::cout << "Unknown topology optimizer, aborting" << std::endl;
			abort();
		}
		return toOptimizer;
	}
}
}

