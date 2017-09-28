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

#include "topoptuniverse.h"
#include "mpihandler.h"
#include "outputhandler.h"

// TO Representation
#include "topoptrep.h"
#include "tomesh.h"
#include "torfactory.h"

// TO Optimizer
#include "topopt.h"
#include "topoptfactory.h"

// TO Objective Function
#include "topoptobjfun.h"
#include <dlfcn.h> // Unix only

#include <string>
#include <fstream>

namespace Topologies{
TopOptUniverse::TopOptUniverse(int argc, char* argv[])
{
	// Load and parse command line input
	std::string fileName;
	bool useMPI = false;
	if (argc == 1)
		fileName = "./topopt.txt";
	else if(argc == 2)
		fileName = argv[1];
	else if(argc == 3)
	{
		std::string opt = argv[1];
		if(opt == "-mpi")
			useMPI = true;
		fileName = argv[2];
	}
	// Parser main input file
	mainParser.parseFile(fileName);

	// Setup MPI first in case objective function initizlizes MPI
	// If upMPIH is not set (i.e. if USE_MPI is not defined), then other functions ignore it and run serially
	#ifdef USE_MPI
	if(useMPI)
	{
		upMPIH = std::unique_ptr<MPIHandler>(new MPIHandler(argc, argv, mainParser.getNumProcessorsPerOFV()));
		if(upMPIH->getSize() == 1)
			upMPIH.release(); // No MPI!
	}
	#endif

	// Parse representation and objective function
	setupTOR();
	setupTOF();

	// Additional MPI setup
	if(upMPIH)
	{
		upMPIH->setTOOF(toofFunc);
		upMPIH->setTOR(torInterface->clone());
	}

	// Setup output
	for(std::vector<InputLoader::Output>::const_iterator oit = mainParser.outputBegin(); oit != mainParser.outputEnd(); ++oit)
		upOutputVec.push_back(std::unique_ptr<OutputHandler>(new OutputHandler(*oit)));

	// Parse optimizer
	setupTOO();
}

TopOptUniverse::~TopOptUniverse()
{
	destroyTOOF(toofFunc);
	dlclose(toobjfunSO);
}

void TopOptUniverse::setupTOR()
{
	torInterface = TopOptRepFactory::createTopOptRep(mainParser.getRepNodeInfo());
}

void TopOptUniverse::setupTOO()
{
	// Create vector of pointers to OutputHandler objects
	std::vector<OutputHandler*> outVec(upOutputVec.size());
	for(std::size_t k = 0; k < upOutputVec.size(); ++k)
		outVec[k] = upOutputVec[k].get();
	// Generate optimizer
	toOptimizer = TopOptFactory::createTopOpt(mainParser.getOptNodeInfo(), toofFunc, outVec, upMPIH.get());
}

void TopOptUniverse::setupTOF()
{
	// Taken from http://www.tldp.org/HOWTO/html_single/C++-dlopen/
	// First load the shared library defined in tofSharedLibFileName
	toobjfunSO = dlopen(mainParser.getTOFSharedLibName().c_str(), RTLD_LAZY);
	if(!toobjfunSO) 
	{
		std::cerr << "Cannot load library: " << dlerror() << std::endl;;
		abort();
	}
	// reset errors
	dlerror();
	// Load objective function creation function
	create_toof* loadTOOF = (create_toof*) dlsym(toobjfunSO, "create");
	const char* dlsym_error = dlerror();
	if(dlsym_error) 
	{
		std::cerr << "Cannot load symbol create: " << dlsym_error << std::endl;
		abort();
	}
	// Load obj fun destroyer function
	destroyTOOF = (destroy_toof*) dlsym(toobjfunSO, "destroy");
	dlsym_error = dlerror();
	if(dlsym_error) 
	{
		std::cerr << "Cannot load symbol destroy: " << dlsym_error << std::endl;
		abort();
	}
	// Load objective function
	toofFunc = loadTOOF(mainParser.getTOFFileName());
}

void TopOptUniverse::setupInitialGuess(TopOptRep* const initialGuess) const
{
	const InputLoader::InitialGuess& theIG = mainParser.getInitialGuess();
	InitialGuessType theIGT = theIG.getInitialGuessType();
	if(theIGT == igtRandom)
		initialGuess->randomize();
	else if(theIGT == igtConstant)
		initialGuess->initialize(theIG.getConstant());
	else if(theIGT == igtConstantWithNoise)
		initialGuess->initialize(theIG.getConstant(), theIG.getRandRange());
	else if(theIGT == igtFile)
		InputLoader::inputTORData(initialGuess, theIG.getFileName());
	else
	{
		std::cerr << "WARNING: Unknown initial guess type, using default constant" << std::endl;
		initialGuess->initialize(1.);
	}
}

void TopOptUniverse::runProblem()
{
	std::unique_ptr<TopOptRep> initialGuess = torInterface->clone();
	setupInitialGuess(initialGuess.get());
	if(upMPIH)
	{
		if(upMPIH->amIRoot())
		{
			result = toOptimizer->optimize(*initialGuess);
			handleOutput(result.get());
		}
		else
			upMPIH->slaveWaitAndProcessFlag(); // Hold slave nodes here for duration
	}
	else
	{
		result = toOptimizer->optimize(*initialGuess);
		handleOutput(result.get());
	}
}

void TopOptUniverse::handleOutput(const TopOptRep* const result) const
{
	for(std::size_t k = 0; k < upOutputVec.size(); ++k)
		upOutputVec[k]->handleOutput(result, toofFunc, true);
}

std::unique_ptr<TopOptRep> TopOptUniverse::getResult() 
{
	return result->clone();
}
}

