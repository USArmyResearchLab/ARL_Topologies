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

#include "topoptrep.h"
#include "topoptchain.h"
#include "topoptfactory.h"
#include "torfactory.h"
#include "mpihandler.h"
#include "inputloaderopt.h"
#include "inputloaderrep.h"
#include "outputwriter.h"
#include "tomesh.h"
#include "csgtree.h"

namespace Topologies{
TopOptChain::TopOptChain(const InputLoader::TOOChain& inputData, TopOptObjFun* inpObjFun, const std::vector<OutputHandler*>& inOutVec, MPIHandler* inpMPIH):
	TopOpt(inpObjFun, inOutVec, inpMPIH),
	oniVec(inputData.oniBegin(), inputData.oniEnd())
{
}

TopOptChain::~TopOptChain()
{
}

std::unique_ptr<TopOptRep> TopOptChain::optimize(const TopOptRep& initialGuess)
{
	std::unique_ptr<TopOptRep> curIG = initialGuess.clone();
	std::unique_ptr<TopOptRep> curRes = initialGuess.clone();
	for(std::size_t k = 0; k < oniVec.size(); ++k)
	{
		if(oniVec[k].getType() == tootRefine)
		{
			std::cout << "Refining TOR" << std::endl;
			refineTOR(curIG.get());
			refineTOR(curRes.get());
		}
		else if(oniVec[k].getType() == tootConvert)
		{
			std::cout << "Converting to LinearSpline" << std::endl;
			curIG = convertTOR(*curIG, oniVec[k]);
			curRes = convertTOR(*curRes, oniVec[k]);
		}
		else if(oniVec[k].getType() == tootContinuation)
		{
			std::cout << "Performing continuation, new TOR tag: " << oniVec[k].getTag() << std::endl;
			curIG = continuationTOR(*curIG, oniVec[k]);
			curRes = continuationTOR(*curRes, oniVec[k]);
		}
		else
		{
			setupTOO(oniVec[k]);
			curRes = toOptimizer->optimize(*curIG);
			curIG = curRes->clone();
		}
	}
	return curRes;
}

void TopOptChain::refineTOR(TopOptRep* torToRefine)
{
	torToRefine->refine();
	if(pMPIH)
	{
		// Must update MPIHandler's TOR so that slaves properly evaluate the new TOR
		std::unique_ptr<TopOptRep> mpiTOR = torToRefine->clone();
		pMPIH->rootConvertTOR(std::move(mpiTOR));
	}
}

std::unique_ptr<TopOptRep> TopOptChain::convertTOR(const TopOptRep& torToConvert, const InputLoader::OptNodeInfo& oni)
{
	// TODO: Currently, this converts to CSGTreeRep, meant for post-processing
	// Fix to be more flexible and able to convert to anything
	// Load input file
	assert(torToConvert.getDimension() == 2);
	InputLoader::RepNodeInfo rni(torToConvert.getDimension() == 2 ? tortCSG2D : tortCSG3D, oni);
	InputLoader::TORCSGTree torParser(rni.getTypeName());
	torParser.parseNode(rni);
	// Get segments
	std::vector<Mesh_Segment_2> segVec;
	torToConvert.get2DSegments(segVec);
	std::unique_ptr<TopOptRep> outTOR = std::unique_ptr<TopOptRep>(new CSGTreeRep(torParser, segVec));
	if(pMPIH)
	{
		// Must update MPIHandler's TOR so that slaves properly evaluate the new TOR
		std::unique_ptr<TopOptRep> mpiTOR = outTOR->clone();
		pMPIH->rootConvertTOR(std::move(mpiTOR), tortCSG2D);
	}
	std::cout << "Converting to CSGTree" << std::endl;
	std::cout << "New CSGTree vol. frac.: " << outTOR->computeVolumeFraction() << std::endl;
	std::unique_ptr<TOMesh> tmpMesh = outTOR->getOutputMesh();
	OutputWriter::plotMeshMatlab(tmpMesh.get(), "csgMesh.m");
	return outTOR;
}

std::unique_ptr<TopOptRep> TopOptChain::continuationTOR(const TopOptRep& torToConvert, const InputLoader::OptNodeInfo& oni)
{
	// Load input params for new TOR
	TORType myType = torToConvert.getType();
	std::unique_ptr<TopOptRep> outTOR = TopOptRepFactory::createTopOptRep(myType, oni);
	// Set new TOR vals to previous TOR
	std::vector<std::vector<int>> discreteVars;
	std::vector<std::vector<double>> realVars;
	torToConvert.getMPIRep(discreteVars, realVars);
	outTOR->setMPIRep(discreteVars, realVars);
	// Update slave nodes' copy of current TOR
	if(pMPIH)
  {
		// Must update MPIHandler's TOR so that slaves properly evaluate the new TOR
		std::unique_ptr<TopOptRep> mpiTOR = outTOR->clone();
		pMPIH->rootConvertTOR(std::move(mpiTOR));
  }
	return outTOR;
}

void TopOptChain::setupTOO(const InputLoader::OptNodeInfo& oni)
{
	toOptimizer = TopOptFactory::createTopOpt(oni, pObjFun, outputVec, pMPIH);
}
}

