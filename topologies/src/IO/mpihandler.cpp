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

#include <iostream>
#include "mpihandler.h"
#include "topoptobjfun.h"
#include "torfactory.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

#define FLAG_DONE 0
#define FLAG_FORM_GROUPS 1
#define FLAG_CLEAR_GROUPS 2
#define FLAG_CONVERT_TOR 3
#define FLAG_FE_ST 4 // Single threaded
#define FLAG_CE_ST 5
#define FLAG_GE_ST 6
#define FLAG_GCE_ST 7
#define FLAG_FANDG_ST 8
#define FLAG_FE_MT 9 // Multi-threaded
#define FLAG_CE_MT 10
#define FLAG_GE_MT 11
#define FLAG_GCE_MT 12
#define FLAG_FANDG_MT 13
#define FLAG_PRINT_MT 14 // Print results using multi-threaded objective func

namespace Topologies{
#ifdef USE_MPI
MPIHandler::MPIHandler(int argc, char* argv[], unsigned inNumProcsPerEval) :	
	myTOOF(nullptr),
	numProcsPerEval(inNumProcsPerEval),
	shouldFinalize(true), // This object initializes MPI, so it should finalize it
	debugPrint(false),
	defaultComm(MPI::COMM_WORLD)
{
	MPI::Init(argc, argv);
	finishSetup();
	std::cout << "MPI initialized on node #" << globalRank << " of " << mpiSize << std::endl;
}

MPIHandler::MPIHandler(MPI::Comm& inComm, unsigned inNumProcsPerEval, const TopOptObjFun* inTOOF, std::unique_ptr<TopOptRep> inTOR) : 
	numProcsPerEval(inNumProcsPerEval),
	shouldFinalize(false), 
	debugPrint(false),
	defaultComm(inComm), 
	myTOOF(inTOOF),
	myTOR(std::move(inTOR))
{
	finishSetup();
}

void MPIHandler::finishSetup()
{
	mpiSize = defaultComm.Get_size();
	globalRank = defaultComm.Get_rank();
	MPI_Barrier(defaultComm);
	localComm = defaultComm;
	communicatorVec.push_back(defaultComm);
	if(debugPrint)
	{
		std::string outFileName;
		std::ostringstream o;
		o << globalRank << ".txt";
		outFileName = "mpiDebugFile" + o.str();
		outFile.open(outFileName.c_str());
		outFile.precision(16);
	}
}

MPIHandler::~MPIHandler()
{
	if(amIRoot())
		rootReleaseSlaves();
	if(debugPrint)
	{
		outFile << "in dtor" << std::endl;
		outFile.close();
	}
	if(shouldFinalize)
	{
		std::cout << "Node #" << globalRank << " calling MPI::Finalize()" << std::endl;
		MPI::Finalize();
	}
}

template<typename Rep>
void MPIHandler::rootBatchEvaluate(const std::vector<Rep>& TORvec, 
	std::vector<std::pair<std::vector<double>, bool>>& resVec, EvalFunction inef)
{
	if(debugPrint) outFile << "numProcsPerEval: " << numProcsPerEval << std::endl;
	if(numProcsPerEval > 1)
	{
		rootBatchEvaluate(TORvec, resVec, numProcsPerEval);
		return;
	}
	MPI::Status status;
	unsigned numGoals = 1;
	unsigned locSize = TORvec.size() < mpiSize ? TORvec.size() : mpiSize;
	for(int kproc = 1; kproc < locSize; kproc++)
	{
		unsigned evalID = kproc - 1;
		rootSetUpForEval(evalID, kproc, inef, true, defaultComm);
		rootSendData(TORvec[evalID], kproc, defaultComm);
	}
	if(debugPrint)
		outFile << "done initial send" << std::endl;
	resVec.resize(TORvec.size());
	for(unsigned keval = locSize; keval <= TORvec.size(); keval++)
	{
		if(debugPrint) outFile << "Calling rootReceiveOneResult" << std::endl;
		// wait for & recieve chromosome information
		unsigned procID = rootReceiveOneResult(resVec);
		// send out new chromosome to be evaluated
		unsigned evalID = keval - 1;
		if(debugPrint) outFile << "Calling rootSetUpForEval with " << evalID << ", " << procID << std::endl;
		rootSetUpForEval(evalID, procID, inef, true, defaultComm);
		if(debugPrint) outFile << "Calling rootSetUpForEval" << std::endl;
		rootSendData(TORvec[evalID], procID, defaultComm);
		if(debugPrint)
			outFile << "Done sending data to process " << procID << std::endl;
	}
	if(debugPrint)
		outFile << "done all sends, waiting for results" << std::endl;
	// receive info from remaining nodes
	for(int kproc = 1; kproc < locSize; kproc++)
		rootReceiveOneResult(resVec);
}

template<typename Rep>
void MPIHandler::rootBatchEvaluate(const std::vector<Rep>& TORvec,
		std::vector<std::pair<std::vector<double>, bool>>& resVec,
		unsigned numProcsPerEval, EvalFunction inef)
{
	// Set up groups
	if(debugPrint) outFile << "Sending flag to form groups" << std::endl;
	// Tell slaves to form groups
	rootSendFlag(FLAG_FORM_GROUPS);
	if(debugPrint) outFile << "Forming groups" << std::endl;
	unsigned requiredSize = numProcsPerEval*TORvec.size() + 1;
	unsigned locSize = requiredSize < mpiSize ? requiredSize : mpiSize; // Number of procs to use
	rootFormGroups(locSize, numProcsPerEval);
	// Begin evaluations
	if(debugPrint) outFile << "Using " << communicatorVec.size() << " groups, and need " << TORvec.size() << " evals" << std::endl;
	std::vector<unsigned> groupIDVec(TORvec.size());
	for(unsigned kgroup = 0; kgroup < communicatorVec.size(); ++kgroup)
	{
		rootSetUpGroupForEval(kgroup, kgroup, inef);
		rootSendDataToGroup(TORvec[kgroup], communicatorVec[kgroup]);
		groupIDVec[kgroup] = kgroup;
	}
	if(debugPrint) outFile << "done initial send" << std::endl;
	// Receive and send
	resVec.resize(TORvec.size());
	for(unsigned keval = communicatorVec.size(); keval < TORvec.size(); ++keval)
	{
		if(debugPrint) outFile << "Calling rootReceiveOneResult" << std::endl;
		// wait for & recieve chromosome information
		unsigned evalID;
		unsigned procID = rootReceiveOneResult(resVec, evalID);
		unsigned kgroup = groupIDVec[evalID];
		if(debugPrint) outFile << "Got evalID: " << evalID << " from proc " << procID << ", in group " << kgroup << std::endl;
		if(debugPrint) outFile << "Calling rootSetUpForEval with " << keval << ", " << kgroup << std::endl;
		rootSetUpGroupForEval(keval, kgroup, inef);
		rootSendDataToGroup(TORvec[keval], communicatorVec[kgroup]);
		// Update groupIDVec
		groupIDVec[keval] = kgroup;
	}
	// Recieve final
	for(unsigned kgroup = 0; kgroup < communicatorVec.size(); ++kgroup)
	{
		if(debugPrint) outFile << "Calling rootReceiveOneResult" << std::endl;
		unsigned evalID;
		unsigned procID = rootReceiveOneResult(resVec, evalID);
		if(debugPrint) outFile << "Got evalID: " << evalID << " from proc " << procID 
														<< ", in group " << groupIDVec[evalID] << std::endl;
	}
	// Disband groups
	if(debugPrint) outFile << "Sending clear groups flag" << std::endl;
	for(auto kComm : communicatorVec)
		rootSendFlagToGroup(FLAG_CLEAR_GROUPS, kComm); // Tell slaves to clear groups
	if(debugPrint) outFile << "Root clearing groups" << std::endl;
	rootClearGroups();
}

void MPIHandler::rootEvaluateDifference(const std::vector<double>& vec1, std::vector<double>& diffVec, std::vector<std::pair<std::vector<double>, bool>>& resVec, EvalFunction inef)
{
	// The 2 vectors vec1 and diffVec represent a realRep and a set of finite differences to compute
	// This avoids having to copy a large number of TopOptRep objects and/or realRep vectors
	MPI::Status status;
	unsigned numGoals = 1;
	unsigned locSize = vec1.size() < mpiSize ? vec1.size() : mpiSize;
	// TODO: Add some check to make sure this matches the MPI rep coming from myTOR
	std::vector<std::vector<double>> realRep;
	std::vector<std::vector<int>> discreteRep;
	realRep.push_back(vec1);
	// Evaluate at x
	rootSetUpForEval(0, 1, inef, true, defaultComm);
	rootSendData(realRep, discreteRep, 1, defaultComm);
	for(int kproc = 2; kproc < locSize; kproc++)
	{
		unsigned evalID = kproc - 1;
		unsigned componentID = evalID - 1;
		rootSetUpForEval(evalID, kproc, inef, true, defaultComm);
		realRep[0][componentID] += diffVec[componentID];
		diffVec[componentID] = realRep[0][componentID] - vec1[componentID];
		rootSendData(realRep, discreteRep, kproc, defaultComm);
		realRep[0][componentID] = vec1[componentID];
	}
	if(debugPrint)
		outFile << "done initial send" << std::endl;
	resVec.resize(vec1.size() + 1);
	for(int keval = locSize; keval <= vec1.size() + 1; keval++)
	{
		if(debugPrint) outFile << "Calling rootReceiveOneResult" << std::endl;
		// wait for & recieve chromosome information
		unsigned procID = rootReceiveOneResult(resVec);
		// send out new chromosome to be evaluated
		unsigned evalID = keval - 1;
		unsigned componentID = evalID - 1;
		if(debugPrint) outFile << "Calling rootSetUpForEval with " << evalID << ", " << procID << std::endl;
		rootSetUpForEval(evalID, procID, inef, true, defaultComm);
		if(debugPrint) outFile << "Calling rootSetUpForEval" << std::endl;
		realRep[0][componentID] += diffVec[componentID];
		diffVec[componentID] = realRep[0][componentID] - vec1[componentID];
		rootSendData(realRep, discreteRep, procID, defaultComm);
		realRep[0][componentID] = vec1[componentID];
		if(debugPrint)
			outFile << "Done sending data to process " << procID << std::endl;
	}
	if(debugPrint)
		outFile << "done all sends, waiting for results" << std::endl;
	// receive info from remaining nodes
	for(int kproc = 1; kproc < locSize; kproc++)
		rootReceiveOneResult(resVec);
}

std::pair<std::vector<double>, bool> MPIHandler::rootEvaluateTOR(const TopOptRep* pTOR, EvalFunction inef)
{
	std::vector<std::pair<std::vector<double>, bool>> resVec(1);
	if(mpiSize == 1 || (numProcsPerEval == 1 && !evalFuncIsGradient(inef))) // Single threaded objective function
		resVec[0] = evaluate(pTOR, inef);
	else // Multi-threaded obj. fun. or gradient
	{
		// Send TOR data
    if(debugPrint) outFile << "Sending info" << std::endl;
		rootSetUpGroupForEval(0, 0, inef); // No groups formed, so communicatorVec will only contain defaultComm
		rootSendDataToGroup(pTOR, communicatorVec.back());
		if(debugPrint)  outFile << "done" << std::endl;
		// In this case, the root will participate in the function evaluation
		// Result is expected only on the root node (since it's the group leader with rank == 0)
		resVec[0] = evaluate(pTOR, inef, communicatorVec.back());
	}
	// Return result
	return resVec[0];
}

void MPIHandler::rootEvaluateTORAndPrint(const TopOptRep* const pTOR, std::string const& fileName)
{
	if(mpiSize == 1 || numProcsPerEval == 1)
		myTOOF->printResult(*pTOR, fileName); // Single threaded
	else
	{
		// Send TOR data
		if(debugPrint) outFile << "Sending info for printResult" << std::endl;
		rootSendFlagToGroup(FLAG_PRINT_MT, communicatorVec.back());
		rootSendDataToGroup(pTOR, communicatorVec.back());
		rootSendStringToGroup(fileName, communicatorVec.back());
		if(debugPrint)  outFile << "done" << std::endl;
		// In this case, the root will participate in the function evaluation
		myTOOF->printResult(*pTOR, fileName, communicatorVec.back());
	}
}

void MPIHandler::rootConvertTOR(std::unique_ptr<TopOptRep> inTOR, TORType inTORT)
{
	if(debugPrint) outFile << "Converting TOR" << std::endl;
	rootSendFlag(FLAG_CONVERT_TOR); // Convert flag
	// Send TOR type
	int buf = (int) inTORT;
	defaultComm.Bcast(&buf, 1, MPI::INT, 0);
	// Send new params
	if(debugPrint) outFile << "Sending new TOR parameters" << std::endl;
	std::vector<std::vector<double>> realParams;
	std::vector<std::vector<int>> discreteParams;
	inTOR->getDefiningParameters(discreteParams, realParams);
	for(int kproc = 1; kproc < mpiSize; kproc++)
		rootSendData(realParams, discreteParams, kproc, defaultComm);
	// Save new TOR
	if(debugPrint) outFile << "Saving new TOR at root node" << std::endl;
	myTOR = std::move(inTOR);
}

void MPIHandler::rootFormGroups(unsigned locSize, unsigned numProcsPerEval)
{
	communicatorVec.clear();
	int numGroups = (locSize - 1)/numProcsPerEval;
	if(debugPrint) outFile << "Setting up " << numGroups << " groups" << std::endl;
	defaultComm.Bcast(&numGroups, 1, MPI_UNSIGNED, 0);
	if(debugPrint) outFile << "Splitting" << std::endl;
	localComm = defaultComm.Split(0, 0); // Assign root to its own group, id 0
	if(debugPrint) outFile << "Done split" << std::endl;
	// Set up intercomms between root and groups and merge into intracomms
	for(Uint k = 1; k <= numGroups; ++k)
	{
		if(debugPrint) outFile << "Setting up Intercomm " << k << std::endl;
		MPI::Intercomm tmpInter = localComm.Create_intercomm(0, defaultComm, k, 0);
		communicatorVec.push_back(tmpInter.Merge(false));
		tmpInter.Free();
		if(debugPrint) outFile << "Done, new communicator size: " << communicatorVec.back().Get_size() << std::endl;
		if(debugPrint) outFile << "Root rank in new comm: " << communicatorVec.back().Get_rank() << std::endl;
	}
}

void MPIHandler::rootClearGroups()
{
	if(debugPrint) outFile << "Clearing all groups" << std::endl;
	// Clear comm vec
	for(auto kComm : communicatorVec)
	{
		if(kComm != defaultComm)
			kComm.Free();
	}
	communicatorVec.clear();
	communicatorVec.push_back(defaultComm);
	// Clear localComm
	if(localComm != defaultComm)
		localComm.Free();
	localComm = defaultComm;
}

void MPIHandler::rootSendString(std::string const& str, unsigned procID, MPI::Comm& communicator) const
{
	if(debugPrint) outFile << "Sending string " << str << std::endl;
	// Send size
	unsigned sz = str.length() + 1;
	communicator.Send(&sz, 1, MPI::UNSIGNED, procID, 0);
	// Send string
	communicator.Send(str.c_str(), sz, MPI::CHAR, procID, 0);
}

void MPIHandler::rootSendData(const TopOptRep* sendTOR, unsigned procID, MPI::Comm& communicator) const
{
	std::vector<std::vector<double>> realRep;
	std::vector<std::vector<int>> discreteRep;
	sendTOR->getMPIRep(discreteRep, realRep);
	if(debugPrint)
	{
		outFile << "Sending " << realRep.size() << " real vectors of sizes: ";
		for(unsigned k = 0; k < realRep.size(); ++k)
			outFile << realRep[k].size() << " ";
		outFile << std::endl;
	}
	rootSendData(realRep, discreteRep, procID, communicator);
}

void MPIHandler::rootSendData(const std::vector<double>& realRep, unsigned procID, MPI::Comm& communicator) const
{
	// This uses the same format as rootSendData with both realRep and discreteRep, so the sizes are still sent
	// Send data sizes:
	unsigned vecSizes[2] = {1, 0};
	communicator.Send(vecSizes, 2, MPI::UNSIGNED, procID, 0);
	// Send all real vectors
	// Send array size
	unsigned vsize = realRep.size();
	communicator.Send(&vsize, 1, MPI::UNSIGNED, procID, 0);
	// Send data
	communicator.Send(realRep.data(), realRep.size(), MPI::DOUBLE, procID, 0);
}

void MPIHandler::rootSendData(const std::vector<std::vector<double>>& realRep, const std::vector<std::vector<int>>& discreteRep, 
	unsigned procID, MPI::Comm& communicator) const
{
	// Send data sizes:
	unsigned vecSizes[2] = {realRep.size(), discreteRep.size()};
	communicator.Send(vecSizes, 2, MPI::UNSIGNED, procID, 0);
	// Send all real vectors
	for(unsigned k = 0; k < vecSizes[0]; ++k)
	{
		// Send array size
		unsigned vsize = realRep[k].size();
		communicator.Send(&vsize, 1, MPI::UNSIGNED, procID, 0);
		// Send data
		communicator.Send(realRep[k].data(), realRep[k].size(), MPI::DOUBLE, procID, 0);
	}
	// Send all discrete vectors
	for(unsigned k = 0; k < vecSizes[1]; ++k)
	{
		// Send array size
		unsigned vsize = discreteRep[k].size();
		communicator.Send(&vsize, 1, MPI::UNSIGNED, procID, 0);
		// Send data
		communicator.Send(discreteRep.data(), discreteRep[k].size(), MPI::INT, procID, 0);
	}
}

void MPIHandler::rootSendStringToGroup(std::string const& str, MPI::Comm& groupComm) const
{
	for(unsigned k = 1; k < groupComm.Get_size(); ++k)
		rootSendString(str, k, groupComm);
}

void MPIHandler::rootSendDataToGroup(const TopOptRep* sendTOR, MPI::Comm& groupComm) const
{
	for(unsigned k = 1; k < groupComm.Get_size(); ++k)
		rootSendData(sendTOR, k, groupComm);
}

void MPIHandler::rootSendDataToGroup(const std::vector<double>& realRep, MPI::Comm& groupComm) const
{
	assert(groupID < communicatorVec.size());
	for(unsigned k = 1; k < groupComm.Get_size(); ++k)
		rootSendData(realRep, k, groupComm);
}

void MPIHandler::rootSendDataToGroup(const std::vector<std::vector<double>>& realRep,
	const std::vector<std::vector<int>>& discreteRep, MPI::Comm& groupComm) const
{
	assert(groupID < communicatorVec.size());
	for(unsigned k = 1; k < groupComm.Get_size(); ++k)
		rootSendData(realRep, discreteRep, k, groupComm);
}

void MPIHandler::rootSendFlagToGroup(unsigned flag, MPI::Comm& groupComm) const
{
	for(unsigned k = 1; k < groupComm.Get_size(); ++k)
		groupComm.Send(&flag, 1, MPI::UNSIGNED, k, 0);
}

void MPIHandler::rootSetUpForEval(unsigned evalID, unsigned procID, EvalFunction inef, bool isST, MPI::Comm& communicator) const
{
	if(debugPrint)
		outFile << "sending info to node " << procID << std::endl;
	unsigned sendBuf = convertEFToFlag(inef, isST);
	communicator.Send(&sendBuf, 1, MPI::UNSIGNED, procID, 0); // eval type flag
	sendBuf = evalID;
	communicator.Send(&sendBuf, 1, MPI::UNSIGNED, procID, 0); // evaluation id
}

void MPIHandler::rootSetUpGroupForEval(unsigned evalID, unsigned groupID, EvalFunction inef)
{
	if(debugPrint)
		outFile << "sending info to group " << groupID << std::endl;
	assert(groupID < communicatorVec.size());
	for(unsigned k = 1; k < communicatorVec[groupID].Get_size(); ++k)
		rootSetUpForEval(evalID, k, inef, false, communicatorVec[groupID]);
}

void MPIHandler::rootSendFlag(unsigned flag) 
{
	for(int kproc = 1; kproc < mpiSize; kproc++)
		defaultComm.Send(&flag, 1, MPI::UNSIGNED, kproc, 0);
}

unsigned MPIHandler::rootReceiveOneResult(std::vector<std::pair<std::vector<double>, bool>>& resVec)
{
	unsigned unused;
	return rootReceiveOneResult(resVec, unused);
}

unsigned MPIHandler::rootReceiveOneResult(std::vector<std::pair<std::vector<double>, bool>>& resVec, unsigned& evalID)
{
	MPI::Status status;
	unsigned evalInfo[3];
	if(debugPrint) outFile << "Waiting for result" << std::endl;
	defaultComm.Recv(evalInfo, 3, MPI::UNSIGNED, MPI::ANY_SOURCE, 0, status);
	if(debugPrint)
    outFile << "Received info from node " << status.Get_source() << std::endl;
  if(debugPrint)
    outFile << "Got evalID: " << evalInfo[0] << ", size of resVec: " << resVec.size() << std::endl;
	unsigned numGoals = evalInfo[2];
	std::unique_ptr<Real[]> ofvs(new Real[numGoals]);
	defaultComm.Recv(ofvs.get(), numGoals, MPI::DOUBLE, status.Get_source(), 0, status);
	if(debugPrint)
	{
		outFile << "Result:";
		for(unsigned k = 0; k < numGoals; ++k)
      outFile << " " << ofvs[k];
		outFile << ", valid: " << evalInfo[1] << std::endl;
	}
	std::vector<double> ofvVec(numGoals);
	for(unsigned k = 0; k < numGoals; ++k)
		ofvVec[k] = ofvs[k];
	evalID = evalInfo[0];
	resVec[evalID] = std::make_pair(ofvVec, evalInfo[1]);
	return status.Get_source();
}

void MPIHandler::rootReleaseSlaves()
{
	for(int kproc = 1; kproc < mpiSize; kproc++)
	{
		unsigned sendBuf = 0;
		defaultComm.Send(&sendBuf, 1, MPI::UNSIGNED, kproc, 0); // eval type flag
	}
}

void MPIHandler::slaveEvaluateAndReturnResult(unsigned evalID, EvalFunction inef, bool isST)
{
	if(!isRunnable())
	{
		std::cerr << "In MPIHandler::slaveEvaluateAndReturnResult, myTOOF or myTOR is uninitialized!" << std::endl;
		abort();
	}
	if(debugPrint)
		outFile << "Evaluating" << std::endl;
	// First evaluate
	std::pair<std::vector<double>, bool> result;
	if(isST)
		result = evaluate(inef);
	else
		result = evaluate(inef, localComm);
	// Send back result
	if(debugPrint)
	{
		outFile << "Done, sending #" << evalID;
		outFile << " result:";
		for(unsigned k = 0; k < result.first.size(); ++k)
			outFile << " " << result.first[k];
		outFile << ", which is " << (result.second ? "" : "not ") << "valid." << std::endl;
	}
	if(isST || (!isST && localComm.Get_rank() == 0)) // Only return result for single threaded or group leader
	{
		unsigned rootID = 0;
		unsigned numGoals = result.first.size(); 
		unsigned evalInfo[3];
		evalInfo[0] = evalID;
		evalInfo[1] = result.second;
		evalInfo[2] = numGoals;
		defaultComm.Send(evalInfo, 3, MPI::UNSIGNED, rootID, 0);
		std::unique_ptr<Real[]> ofvs(new Real[numGoals]);
		for(unsigned k = 0; k < numGoals; ++k)
			ofvs[k] = result.first[k];
		defaultComm.Send(ofvs.get(), numGoals, MPI::DOUBLE, rootID, 0);
	}
}

void MPIHandler::slaveFormGroups()
{
	// Split into a given number of groups
	if(debugPrint) outFile << "Setting up groups" << std::endl;
	unsigned numGroups;
	defaultComm.Bcast(&numGroups, 1, MPI_UNSIGNED, 0);
	unsigned groupID = defaultComm.Get_rank() % numGroups + 1;
	if(debugPrint) outFile << "Got num groups: " << numGroups << ", my id: " << groupID << std::endl;
	localComm = defaultComm.Split(groupID, 0);
	if(debugPrint) outFile << "Done split" << std::endl;
	// Set up intercomms to communicate with root node.
	MPI::Intercomm tmpInter = localComm.Create_intercomm(0, defaultComm, 0, 0);
	communicatorVec.push_back(tmpInter.Merge(true));
	if(debugPrint) outFile << "Done, local size: " << localComm.Get_size() << ", rank: " << localComm.Get_rank() << std::endl;
	if(debugPrint) outFile << "Done, inter size: " << communicatorVec.back().Get_size() << ", rank: " << communicatorVec.back().Get_rank() << std::endl;
}

void MPIHandler::slaveClearGroups()
{
	rootClearGroups(); // Same procedure
}

std::string MPIHandler::slaveReceiveString()
{
	MPI::Intracomm comm = communicatorVec.empty() ? defaultComm : communicatorVec.back();
	// Receive size
	unsigned sz;
	MPI::Status status;
	comm.Recv(&sz, 1, MPI::UNSIGNED, MPI::ANY_SOURCE, 0, status);
	// Receive string
	std::unique_ptr<char[]> cStr(new char[sz]);
	comm.Recv(cStr.get(), sz, MPI::CHAR, MPI::ANY_SOURCE, 0, status);
	// Return string
	return std::string(cStr.get());
}

void MPIHandler::checkRunnable() const
{
	if(!isRunnable())
	{
		std::cerr << "In MPIHandler::slaveEvaluateAndReturnResult, myTOOF or myTOR is uninitialized!" << std::endl;
		abort();
	}
}

unsigned MPIHandler::slaveReceiveEvalID()
{
	if(debugPrint) outFile << "Waiting for TOR eval id" << std::endl;
	// Evaluation ID so that root knows which TOR is being evaluated on return
	MPI::Status status;
	unsigned evalID;
	MPI::Intracomm comm = communicatorVec.empty() ? defaultComm : communicatorVec.back();
	comm.Recv(&evalID, 1, MPI::UNSIGNED, MPI::ANY_SOURCE, 0, status);
	return evalID;
}

void MPIHandler::slaveReceiveDataAndSetupTOR()
{
	if(debugPrint) outFile << "Waiting for TOR data" << std::endl;
	// Recv data
	std::vector<std::vector<double>> realRep;
	std::vector<std::vector<int>> discreteRep;
	slaveReceiveData(realRep, discreteRep);
	// Finish set up
	if(debugPrint) outFile << "Got TOR, setting up" << std::endl;
	myTOR->setMPIRep(discreteRep, realRep);
}

void MPIHandler::slaveReceiveData(std::vector<std::vector<double>>& realRep, std::vector<std::vector<int>>& discreteRep)
{
	// Receive data sizes:
	MPI::Status status;
	unsigned vecSizes[2];
	communicatorVec.back().Recv(vecSizes, 2, MPI::UNSIGNED, MPI::ANY_SOURCE, 0, status);
	// Receive all real vectors
	realRep.resize(vecSizes[0]);
	for(unsigned k = 0; k < vecSizes[0]; ++k)
	{
		// Receive array size
		unsigned vsize;
		communicatorVec.back().Recv(&vsize, 1, MPI::UNSIGNED, MPI::ANY_SOURCE, 0, status);
		// Receive data
		realRep[k].resize(vsize);
		communicatorVec.back().Recv(realRep[k].data(), vsize, MPI::DOUBLE, MPI::ANY_SOURCE, 0, status);
		if(debugPrint)
		{
			outFile << "Received data: ";
			for(auto val : realRep[k])
				outFile << val << " ";
			outFile << std::endl;
		}
	}
	// Recv all discrete vectors
	discreteRep.resize(vecSizes[1]);
	for(unsigned k = 0; k < vecSizes[1]; ++k)
	{
		// Receive array size
		unsigned vsize;
		communicatorVec.back().Recv(&vsize, 1, MPI::UNSIGNED, MPI::ANY_SOURCE, 0, status);
		// Receive data
		discreteRep[k].resize(vsize);
		communicatorVec.back().Recv(discreteRep[k].data(), vsize, MPI::INT, MPI::ANY_SOURCE, 0, status);
	}
}

void MPIHandler::slaveReceiveNewTORParams()
{
	if(debugPrint) outFile << "Converting TOR" << std::endl;
  // Recv TOR type	
	int tortype;
	defaultComm.Bcast(&tortype, 1, MPI::INT, 0);
  // Recv new params
	std::vector<std::vector<double>> realParams;
	std::vector<std::vector<int>> discreteParams;
	slaveReceiveData(realParams, discreteParams);
	// Set up new TOR
	setupNewTOR((TORType)tortype, realParams, discreteParams);
}

void MPIHandler::slavePrintResult()
{
	checkRunnable();
	slaveReceiveDataAndSetupTOR();
	std::string fileName = slaveReceiveString();
	myTOOF->printResult(*myTOR, fileName, localComm);
}

void MPIHandler::slaveEvaluate(unsigned flag)
{
	checkRunnable();
	unsigned evalID = slaveReceiveEvalID();
	slaveReceiveDataAndSetupTOR();
	EvalFunction curef = convertFlagToEF(flag);
	slaveEvaluateAndReturnResult(evalID, curef, isEvalSingleThreaded(flag));
}

void MPIHandler::slaveWaitAndProcessFlag()
{
	unsigned flag;
	MPI::Status status;
	if(debugPrint) outFile << "Waiting for flag, current group size: " << communicatorVec.back().Get_size() << ", rank: " << communicatorVec.back().Get_rank() << std::endl;
	communicatorVec.back().Recv(&flag, 1, MPI::UNSIGNED, MPI::ANY_SOURCE, 0, status);
	if(debugPrint) outFile << "Got flag " << flag << std::endl;
	while(flag != FLAG_DONE)
	{
		if(flagIsEvaluation(flag)) // Evaluate objective function
			slaveEvaluate(flag);
		else if(flag == FLAG_FORM_GROUPS) // Set up for multi-threaded OF
			slaveFormGroups();
		else if(flag == FLAG_CLEAR_GROUPS)
			slaveClearGroups();
		else if(flag == FLAG_CONVERT_TOR)
			slaveReceiveNewTORParams();
		else if(flag == FLAG_PRINT_MT)
			slavePrintResult();
		if(debugPrint) outFile << "Waiting for flag, current group size: " << communicatorVec.back().Get_size() << ", rank: " << communicatorVec.back().Get_rank() << std::endl;
		communicatorVec.back().Recv(&flag, 1, MPI::UNSIGNED, MPI::ANY_SOURCE, 0, status);
		if(debugPrint) outFile << "Got flag " << flag << std::endl;
	}
}

void MPIHandler::setupNewTOR(TORType tortype, const std::vector<std::vector<double>>& realRep,
                  						const std::vector<std::vector<int>>& discreteRep)
{
	myTOR = TopOptRepFactory::createTopOptRep(tortype, realRep, discreteRep);
}

bool MPIHandler::isRunnable() const
{
	return (myTOOF != nullptr) && myTOR;
}

EvalFunction MPIHandler::convertFlagToEF(unsigned flag) const
{
	if(flag == FLAG_FE_ST || flag == FLAG_FE_MT)
		return efF;
	else if(flag == FLAG_CE_ST || flag == FLAG_CE_MT)
		return efC;
	else if(flag == FLAG_GE_ST || flag == FLAG_GE_MT)
		return efG;
	else if(flag == FLAG_GCE_ST || flag == FLAG_GCE_MT)
		return efGC;
	else if(flag == FLAG_FANDG_ST || flag == FLAG_FANDG_MT)
		return efFandG;
	return efF; // Default case
}

unsigned MPIHandler::convertEFToFlag(EvalFunction inEF, bool isST) const
{
	if(inEF == efF)
		return isST ? FLAG_FE_ST : FLAG_FE_MT;
	else if(inEF == efC)
		return isST ? FLAG_CE_ST : FLAG_CE_MT;
	else if(inEF == efG)
		return isST ? FLAG_GE_ST : FLAG_GE_MT;
	else if(inEF == efGC)
		return isST ? FLAG_GCE_ST : FLAG_GCE_MT;
	else if(inEF == efFandG)
		return isST ? FLAG_FANDG_ST : FLAG_FANDG_MT;
	return FLAG_FE_ST; // default case
}

bool MPIHandler::evalFuncIsGradient(EvalFunction inEF) const
{
	return inEF == efG || inEF == efGC;
}

bool MPIHandler::flagIsEvaluation(unsigned flag) const
{
	return flag == FLAG_FE_ST || flag == FLAG_CE_ST || flag == FLAG_GE_ST || flag == FLAG_GCE_ST || flag == FLAG_FANDG_ST ||
					flag ==	FLAG_FE_MT || flag == FLAG_CE_MT || flag == FLAG_GE_MT || flag == FLAG_GCE_MT || flag == FLAG_FANDG_MT;
}

bool MPIHandler::isEvalSingleThreaded(unsigned flag) const
{
	return flag == FLAG_FE_ST || flag == FLAG_CE_ST || flag == FLAG_GE_ST || flag == FLAG_GCE_ST || flag == FLAG_FANDG_ST;
}

#else
// Non-mpi versions
// These exist to reduce the occurance of ifdef USE_MPI throughout the code
// Some functions do nothing, but some implement serial versions of the above.
MPIHandler::MPIHandler(int argc, char* argv[], unsigned inNumProcsPerEval) :
	myTOOF(nullptr),
  numProcsPerEval(inNumProcsPerEval),
  debugPrint(false),
	shouldFinalize(true),
  mpiSize(1),
  globalRank(0)
{
}

MPIHandler::~MPIHandler()
{
}

void MPIHandler::rootReleaseSlaves()
{
}

void MPIHandler::slaveWaitAndProcessFlag()
{
	// This should never be called in non-MPI mode
}

template<typename Rep>
void MPIHandler::rootBatchEvaluate(const std::vector<Rep>& TORvec, 
																	std::vector<std::pair<std::vector<double>, bool>>& resVec, 
																	EvalFunction inef)
{
	resVec.resize(TORvec.size());
	for(unsigned k = 0; k < TORvec.size(); ++k)
		resVec[k] = evaluate(TORVec[k], inef);
}

template<typename Rep>
void MPIHandler::rootBatchEvaluate(const std::vector<Rep>& TORvec, 
																	std::vector<std::pair<std::vector<double>, bool>>& resVec, 
																	unsigned numProcsPerEval, EvalFunction inef)
{
	rootBatchEvaluate(TORvec, resVec, inef);
}

std::pair<std::vector<double>, bool> MPIHandler::rootEvaluateTOR(const TopOptRep* pTOR, EvalFunction inef)
{
	std::pair<std::vector<double>, bool> res;
	rootBatchEvaluate({pTOR}, resVec, inef);
	return res;
}

void MPIHandler::rootEvaluateTORAndPrint(const TopOptRep* const pTOR, std::string const& fileName)
{
	pTOR->printResult(fileName);
}

void MPIHandler::rootEvaluateDifference(const std::vector<double>& vec1, std::vector<double>& diffVec, 
																				std::vector<std::pair<std::vector<double>, bool>>& resVec, EvalFunction inef)
{
  // The 2 vectors vec1 and diffVec represent a realRep and a set of finite differences to compute
  // This avoids having to copy a large number of TopOptRep objects and/or realRep vectors
  unsigned numGoals = 1;
  unsigned locSize = vec1.size() < mpiSize ? vec1.size() : mpiSize;
  std::vector<std::vector<double>> realRep;
  std::vector<std::vector<int>> discreteRep;
  realRep.push_back(vec1);
	resVec.resize(vec1.size() + 1);
	// Evaluate at x
	resVec[0] = evaluate(realRep, discreteRep, inef);
	for(unsigned k = 0; k < vec1.size(); ++k)
	{
		// Set up difference in realRep
		realRep[0][k] += diffVec[k];
		diffVec[k] = realRep[0][k] - vec1[k];
		// Evaluate
		resVec[k+1] = evaluate(realRep, discreteRep, inef);
		// Reset realRep
		realRep[0][k] = vec1[k];
	}
}

void MPIHandler::rootConvertTOR(std::unique_ptr<TopOptRep> inTOR, TORType inTORT)
{
	myTOR = std::move(inTOR);
}

#endif

std::pair<std::vector<double>, bool> MPIHandler::evaluate(std::vector<std::vector<double>>& realRep,
		std::vector<std::vector<int>>& discreteRep, EvalFunction inef)
{
	myTOR->setMPIRep(discreteRep, realRep);
	return evaluate(inef);
}

std::pair<std::vector<double>, bool> MPIHandler::evaluate(std::vector<double> const& realRep, EvalFunction inef)
{
	myTOR->setMPIRep({}, {realRep});
	return evaluate(inef);
}

std::pair<std::vector<double>, bool> MPIHandler::evaluate(EvalFunction inef)
{
	return evaluate(myTOR.get(), inef);
}

std::pair<std::vector<double>, bool> MPIHandler::evaluate(const TopOptRep* pTOR, EvalFunction inef)
{
	std::pair<std::vector<double>, bool> result;
	if(inef == efF)
		myTOOF->f(*pTOR, result);
	else if(inef == efC)
		myTOOF->c(*pTOR, result);
	else if(inef == efG)
		myTOOF->g(*pTOR, result);
	else if(inef == efFandG)
	{
		std::pair<std::vector<double>, bool> gRes;
		myTOOF->fAndG(*pTOR, result, gRes);
		// Append the vectors and set error flag
		result.first.reserve(result.first.size() + gRes.first.size());
		result.first.insert(result.first.end(), gRes.first.begin(), gRes.first.end());
		result.second &= gRes.second;
	}
	else // gc
		myTOOF->gc(*pTOR, result);
	return result;
}

std::pair<std::vector<double>, bool> MPIHandler::evaluate(EvalFunction inef, MPI::Comm& inComm)
{
	return evaluate(myTOR.get(), inef, inComm);
}

std::pair<std::vector<double>, bool> MPIHandler::evaluate(const TopOptRep* pTOR, EvalFunction inef, MPI::Comm& inComm)
{
	std::pair<std::vector<double>, bool> result;
	if(inef == efF)
		myTOOF->f(*pTOR, result, inComm);
	else if(inef == efC)
		myTOOF->c(*pTOR, result, inComm);
	else if(inef == efG)
		myTOOF->g(*pTOR, result, inComm);
	else if(inef == efFandG)
	{
		std::pair<std::vector<double>, bool> gRes;
		myTOOF->fAndG(*pTOR, result, gRes, inComm);
		// Append the vectors and set error flag
		result.first.reserve(result.first.size() + gRes.first.size());
		result.first.insert(result.first.end(), gRes.first.begin(), gRes.first.end());
		result.second &= gRes.second;
	}
	else // efGC
		myTOOF->gc(*pTOR, result, inComm);
	return result;
}

void MPIHandler::rootConvertTOR(std::unique_ptr<TopOptRep> inTOR)
{
	TORType tmp = inTOR->getType();
	rootConvertTOR(std::move(inTOR), tmp);
}

void MPIHandler::copyFirstObjective(std::vector<std::pair<double, bool>>& resVec,
	const std::vector<std::pair<std::vector<double>, bool>>& resVec2) const
{
	resVec.resize(resVec2.size());
	for(unsigned k = 0; k < resVec.size(); ++k)
	{
		double ofv = 0.;
		if(!resVec2[k].first.empty())
			ofv = resVec2[k].first[0];
		resVec[k] = std::make_pair(ofv, resVec2[k].second);
	}
}

template<typename Rep>
void MPIHandler::rootBatchEvaluate(const std::vector<Rep>& TORvec,
                                    std::vector<std::pair<double, bool>>& resVec, EvalFunction inef)
{
	std::vector<std::pair<std::vector<double>, bool>> tmpResVec;
	rootBatchEvaluate(TORvec, tmpResVec, inef);
	copyFirstObjective(resVec, tmpResVec);
}

template<typename Rep>
void MPIHandler::rootBatchEvaluate(const std::vector<Rep>& TORvec,
                                  std::vector<std::pair<double, bool>>& resVec, unsigned numProcsPerEval,
                                  EvalFunction inef)
{
	std::vector<std::pair<std::vector<double>, bool>> tmpResVec;
	rootBatchEvaluate(TORvec, tmpResVec, numProcsPerEval, inef);
	copyFirstObjective(resVec, tmpResVec);
}

void MPIHandler::rootEvaluateDifference(const std::vector<double>& vec1, std::vector<double>& diffVec,
                                    std::vector<std::pair<double, bool>>& resVec, EvalFunction inef)
{
	std::vector<std::pair<std::vector<double>, bool>> tmpResVec;
	rootEvaluateDifference(vec1, diffVec, tmpResVec, inef);
	copyFirstObjective(resVec, tmpResVec);
}

// Explicit instantiations
template void MPIHandler::rootBatchEvaluate(const std::vector<TopOptRep*>&, std::vector<std::pair<double, bool>>&, EvalFunction);
template void MPIHandler::rootBatchEvaluate(const std::vector<TopOptRep*>&, std::vector<std::pair<double, bool>>&, 
	unsigned, EvalFunction);
template void MPIHandler::rootBatchEvaluate(const std::vector<TopOptRep*>&, std::vector<std::pair<std::vector<double>, bool>>&, 
	EvalFunction);
template void MPIHandler::rootBatchEvaluate(const std::vector<TopOptRep*>&, std::vector<std::pair<std::vector<double>, bool>>&, 
	unsigned, EvalFunction);
template void MPIHandler::rootBatchEvaluate(const std::vector<std::vector<double>>&, std::vector<std::pair<double, bool>>&, 
	EvalFunction);
template void MPIHandler::rootBatchEvaluate(const std::vector<std::vector<double>>&, std::vector<std::pair<double, bool>>&, 
	unsigned, EvalFunction);
template void MPIHandler::rootBatchEvaluate(const std::vector<std::vector<double>>&, 
	std::vector<std::pair<std::vector<double>, bool>>&,	EvalFunction);
template void MPIHandler::rootBatchEvaluate(const std::vector<std::vector<double>>&, 
	std::vector<std::pair<std::vector<double>, bool>>&, unsigned, EvalFunction);
}


