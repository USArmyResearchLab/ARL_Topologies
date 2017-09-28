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

#ifndef TOPOPTUNIVERSE_H
#define TOPOPTUNIVERSE_H

// Class to manage optimization process.
// Responsible for reading input files and creating appropriate TOO, TOR, and TOOF objects

#include <memory>
#include "topoptobjfun.h"
#include "inputloadertopopt.h"

namespace Topologies{
class TopOptRep;
class TopOpt;
class MPIHandler;
class Parser;
class OutputHandler;

//! The main class for Topologies, TopOptUniverse handles input, output, and running of the optimizers
/*! This is the main class for Topologies and is responsible for parsing input, creating all 
 *  necessary objects, setting up the output, setting up MPI, and running the optimizers.  
 *  It's public interface is rather small, and only consists of a constructor whihc takes the
 *  command line arguments (which should contain the intput file to use) and a runProblem function.
 */
class TopOptUniverse
{
public:
	//! Constructor that takes the command line arguments
	TopOptUniverse(int argc, char* argv[]);
	~TopOptUniverse();

	//! Runs the optimization problem
	void runProblem();
	//! Returns a unique_ptr to the result
	std::unique_ptr<TopOptRep> getResult();
private:
	void setupTOR();
	void setupTOO();
	void setupTOF();
	void setupInitialGuess(TopOptRep* const initialGuess) const;
	void handleOutput(const TopOptRep* const result) const;

	std::unique_ptr<TopOpt> toOptimizer;
	std::unique_ptr<TopOptRep> torInterface, result;
	TopOptObjFun* toofFunc; // Memory is managed in a shared library
	std::unique_ptr<MPIHandler> upMPIH;
	std::vector<std::unique_ptr<OutputHandler>> upOutputVec;
	InputLoader::TopologiesParser mainParser;
	// Shared library
	void* toobjfunSO;
	destroy_toof* destroyTOOF;
private:
	TopOptUniverse& operator=(TopOptUniverse);
	TopOptUniverse(const TopOptUniverse&);
	TopOptUniverse(TopOptUniverse&&);
};
}
#endif

