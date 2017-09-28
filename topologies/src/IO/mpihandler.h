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

#ifndef MPIHANDLER_H
#define MPIHANDLER_H

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <vector>
#include <fstream>
#include <memory>
#include "cgal_types.h"
#include "topologiesdefs.h"
#include "topoptrep.h"

namespace Topologies{
class TopOptObjFun;

//! A namespace for converting STL containers to C-style arrays for MPI usage
/*! Note that these will return arrays with a size of N+1, if the input vector is size N.
*  The last element is a check sum, sum of the absolute values of the entries.
*/
namespace MPIHelperFunctions
{
	std::unique_ptr<Real[]> packData(const std::vector<Real>& realVec);
	std::unique_ptr<int[]> packData(const std::vector<int>& intVec);
	void unpackData(const int *const inArray, unsigned arraySize, std::vector<int>& intVec);
	void unpackData(const Real *const inArray, unsigned arraySize, std::vector<Real>& realVec);
	bool checkArrayCheckSum(const int *const inArray, unsigned arraySize);
	bool checkArrayCheckSum(const Real *const inArray, unsigned arraySize);
}

//! A wrapper class for MPI parallization
/*! This class is responsible for implementing all parallelization at the optimizer level.
*  Implementations of TopOptObjFun can also be parallelized, but they are responsible
*  for their own parallel implemenation through a communicator passed by an MPIHandler object.
*  A preprocessor definition USE_MPI can be turned off to remove all MPI code and make
*  compilation in a non-MPI environment possible.  In that case, this class will do nothing, but
*  return 1 for its MPI size and 0 for its MPI rank.
*/
class MPIHandler
{
public:
	//! Constructor which takes the command line arguments and the number of processors to use for objective function evaluations
	/*! This is the only constructor as only one MPIHandler object should be in existance at a time.  This constructor calls
	*  MPI::Init and saves this process's rank and total number of processors.
	*/
	MPIHandler(int argc, char* argv[], unsigned inNumProcsPerEval);
	MPIHandler(MPI::Comm& inComm, unsigned inNumProcsPerEval, const TopOptObjFun* inTOOF, std::unique_ptr<TopOptRep> inTOR);
	~MPIHandler();

	//! Release slave nodes from servitude
	/*! For MPI parallelized runs, all nodes but 1 enter a loop where they receive flags from the root process
	 *  and perform evaluations.  This will release them so that the program can exit.
	 */
	void rootReleaseSlaves();
	//! Main slave node function where the slave nodes receive flags from the root node
	void slaveWaitAndProcessFlag();
	//! Returns the current TopOptObjFun
	const TopOptObjFun* getTOOF() const {return myTOOF;}
	//! Sets the current TopOptObjFun to that specified in `inTOOF`
	void setTOOF(const TopOptObjFun* inTOOF) {myTOOF = inTOOF;}
	//! Sets the current TopOptRep to that in `inTOR`.  The MPIHandler object assumes control of the `inTOR`.
	void setTOR(std::unique_ptr<TopOptRep> inTOR) {myTOR = std::move(inTOR);}
	//! Direct slave nodes to evalute a group of TopOptRep objects
	void rootBatchEvaluate(const std::vector<TopOptRep*>& TORvec, std::vector<std::pair<double, bool>>& resVec, EvalFunction inef = efF);
	//! Direct slave nodes to evalute a group of TopOptRep objects, using parallel objective function evaluations
	void rootBatchEvaluate(const std::vector<TopOptRep*>& TORvec, std::vector<std::pair<double, bool>>& resVec, unsigned numProcsPerEval, EvalFunction inef = efF);
	//! Direct slave nodes to evalute a group of TopOptRep objects
	void rootBatchEvaluate(const std::vector<TopOptRep*>& TORvec, std::vector<std::pair<std::vector<double>, bool>>& resVec, EvalFunction inef = efF);
	//! Direct slave nodes to evalute a group of TopOptRep objects, using parallel objective function evaluations
	void rootBatchEvaluate(const std::vector<TopOptRep*>& TORvec, std::vector<std::pair<std::vector<double>, bool>>& resVec, unsigned numProcsPerEval, EvalFunction inef = efF);
	//! Compute, in parallel, the finite difference approximation of the gradient of the current TopOptRep.
	/*! This function specifies the continuous representation of a TopOptRep in vec1, and uses the values in diffVec to comptue
	 *  a finite difference approximation of the gradient of the objective function with respect to each parameter in vec1
	 */
	void rootEvaluateDifference(const std::vector<double>& vec1, std::vector<double>& diffVec,
		std::vector<std::pair<std::vector<double>, bool>>& resVec, EvalFunction inef = efF);
	//! Compute, in parallel, the finite difference approximation of the gradient of the current TopOptRep.
	/*! This function specifies the continuous representation of a TopOptRep in vec1, and uses the values in diffVec to comptue
	 *  a finite difference approximation of the gradient of the objective function with respect to each parameter in vec1
	 */
	void rootEvaluateDifference(const std::vector<double>& vec1, std::vector<double>& diffVec,
		std::vector<std::pair<double, bool>>& resVec, EvalFunction inef = efF);
	//! Evaluates one TOR using a parallelized objective function
	std::pair<std::vector<double>, bool> rootEvaluateTOR(const TopOptRep* pTOR, EvalFunction inef = efF);
	//! Converts the current TopOptRep to that specified by `torName` and `inTOR`
	void rootConvertTOR(std::unique_ptr<TopOptRep> inTOR, TORType inTORT);
	//! Converts the current TopOptRep to that specified by `torName` and `inTOR`
	void rootConvertTOR(std::unique_ptr<TopOptRep> inTOR);
	//! Returns whether or not this is the root node
	bool amIRoot() const;
	//! Returns the MPI rank
	unsigned getRank() const;
	//! REturns the number of MPI nodes
	unsigned getSize() const;

private:
#ifdef USE_MPI
	void finishSetup();
	bool isRunnable() const;
	void rootFormGroups(unsigned locSize, unsigned numProcsPerEval);
	void rootClearGroups();
	void rootSendData(const TopOptRep* sendTOR, unsigned procID, MPI::Comm& communicator);
	void rootSendData(const std::vector<std::vector<double>>& realRep, 
		const std::vector<std::vector<int>>& discreteRep, unsigned procID, MPI::Comm& communicator);
	void rootSendDataToGroup(const TopOptRep* sendTOR, MPI::Comm& groupComm);
	void rootSendDataToGroup(const std::vector<std::vector<double>>& realRep,
		const std::vector<std::vector<int>>& discreteRep, MPI::Comm& groupComm);
	void rootSetUpForEval(unsigned evalID, unsigned procID, EvalFunction inef, bool isST, MPI::Comm& communicator);
	void rootSetUpGroupForEval(unsigned evalID, unsigned groupID, EvalFunction inef);
	unsigned rootReceiveOneResult(std::vector<std::pair<std::vector<double>, bool>>& resVec);
	unsigned rootReceiveOneResult(std::vector<std::pair<std::vector<double>, bool>>& resVec, unsigned& evalID);
	void rootSendFlag(unsigned flag) const;
	void rootSendFlagToGroup(unsigned flag, MPI::Comm& groupComm) const;

	void slaveFormGroups();
	void slaveClearGroups();
	unsigned slaveReceiveDataAndSetupTOR();
	void slaveReceiveData(std::vector<std::vector<double>>& realRep, std::vector<std::vector<int>>& discreteRep);
	void slaveEvaluateAndReturnResult(unsigned evalID, EvalFunction inef, bool isST);
	void slaveReceiveNewTORParams();

	void setupNewTOR(TORType tortype, const std::vector<std::vector<double>>& realRep, 
		const std::vector<std::vector<int>>& discreteRep);

	EvalFunction convertFlagToEF(unsigned flag) const;
	unsigned convertEFToFlag(EvalFunction inEF, bool isST = true) const;
	bool flagIsEvaluation(unsigned flag) const;
	bool evalFuncIsGradient(EvalFunction inEF) const;
	bool isEvalSingleThreaded(unsigned flag) const;

	MPI::Intracomm defaultComm;
	MPI::Intracomm localComm;
	std::vector<MPI::Intracomm> communicatorVec;
#endif
	std::pair<std::vector<double>, bool> evaluate(std::vector<std::vector<double>>& realRep, 
		std::vector<std::vector<int>>& discreteRep, EvalFunction inef);
	std::pair<std::vector<double>, bool> evaluate(EvalFunction inef);
	std::pair<std::vector<double>, bool> evaluate(const TopOptRep* pTOR, EvalFunction inef);
	std::pair<std::vector<double>, bool> evaluate(EvalFunction inef, MPI::Comm& inComm);
	std::pair<std::vector<double>, bool> evaluate(const TopOptRep* pTOR, EvalFunction inef, MPI::Comm& inComm);

	void copyFirstObjective(std::vector<std::pair<double, bool>>& resVec,
		const std::vector<std::pair<std::vector<double>, bool>>& resVec2) const;

	const TopOptObjFun* myTOOF;
	std::unique_ptr<TopOptRep> myTOR;
	std::ofstream outFile;
	bool debugPrint, shouldFinalize;
	unsigned numProcsPerEval;
	unsigned mpiSize, globalRank;
private:
	MPIHandler& operator=(MPIHandler);
	MPIHandler(const MPIHandler&);
	MPIHandler(MPIHandler&&);
};

inline
bool MPIHandler::amIRoot() const
{
	return globalRank == 0;
}

inline
unsigned MPIHandler::getRank() const
{
	return globalRank;
}

inline
unsigned MPIHandler::getSize() const
{
	return mpiSize;
}
}
#endif

