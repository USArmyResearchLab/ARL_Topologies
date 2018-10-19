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

#include "topoptobjfun.h"
#include "topoptrep.h"
#include "mpihandler.h"
#include <memory>
// DEBUG
#include <string>
#include <sstream>
#include "outputwriter.h"
#include "tomesh.h"

namespace Topologies{
void TopOptObjFun::g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	// Default finite difference approximation
	std::vector<double> y;
  inTOR.getRealRep(y);
	std::vector<double> x = y;
	double sqrtEta = sqrt(getEps());
	std::unique_ptr<TopOptRep> locTOR = inTOR.clone();
	std::pair<double, bool> fc = (*this)(*locTOR);
	std::vector<double> outGrad(y.size());
	bool valid = fc.second;
	bool debugPrint = true;
	if(debugPrint)
	{
		std::cout << "In serial gradient computation, need " << y.size() << " components" << std::endl;
		std::cout.precision(16);
		std::cout << "Current ofv: " << fc.first << std::endl;
	}
  for (int j = 0; j < y.size() && valid; ++j)
  {
    double stepSizeJ = sqrtEta*MAX(fabs(x[j]), inTOR.getDataMagnitude()) * (x[j] >= 0 ? 1. : -1.);
		if(x[j] == 1.) // Switch to backward difference, 1 is assumed max for all TORs
			stepSizeJ *= -1.;
    double tempJ = x[j];
    x[j] += stepSizeJ;
    stepSizeJ = x[j] - tempJ;
		locTOR->setRealRep(x.begin(), x.end());
    std::pair<double, bool> fj = (*this)(*locTOR);
		if(debugPrint)
		{
			std::cout << "Done " << j << ", res: " << fj.first << ", dif: " << (fj.first - fc.first) 
								<< ", stepSize: " << stepSizeJ << std::endl;
		}
		outGrad[j] = (fj.first - fc.first)/stepSizeJ;
		valid &= fj.second;
    x[j] = tempJ;
  }
	if(!valid)
		std::cout << "Warning: An objective function computation was flagged as invalid in the gradient" << std::endl;
	outRes = std::make_pair(outGrad, valid);
}

void TopOptObjFun::g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	computeGradient(efF, inTOR, outRes, communicator);
}

void TopOptObjFun::forwardDiff(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPIHandler* pMPIH, EvalFunction inef) const
{
	// default finite difference approximation
	std::vector<double> y;
	inTOR.getRealRep(y);
	std::vector<std::unique_ptr<TopOptRep> > torUPVec;
	torUPVec.reserve(y.size() + 1);
	std::vector<TopOptRep*> torVec(y.size() + 1);
	torUPVec.push_back(inTOR.clone());
	torVec[0] = torUPVec[0].get();
	std::vector<double> stepSizes(y.size());
	double sqrtEta = sqrt(getEps());
	for(std::size_t j = 0; j < y.size() ; ++j)
	{
		double stepSizeJ = sqrtEta*MAX(fabs(y[j]), inTOR.getDataMagnitude()) * (y[j] >= 0 ? 1. : -1.);
		if(y[j] == 1.) // Switch to backward difference, 1 is assumed max for all TORs
			stepSizeJ *= -1.;
		double tempJ = y[j];
		y[j] += stepSizeJ;
		stepSizes[j] = y[j] - tempJ;
		torUPVec.push_back(inTOR.clone());
		torVec[j + 1] = torUPVec.back().get();
		torVec[j + 1]->setRealRep(y.begin(), y.end());
		y[j] = tempJ;
	}
	// MPI evaluate all TORs
	std::vector<std::pair<double, bool> > resVec;
	pMPIH->rootBatchEvaluate(torVec, resVec, inef);
	// Construct gradient from results
	std::vector<double> outGrad(inTOR.getDataSize());
	bool valid = false;
	if(!resVec.empty())
		valid = resVec[0].second;
	for(std::size_t k = 1; k < resVec.size(); ++k)
  {
    valid &= resVec[k].second;
    outGrad[k - 1] = (resVec[k].first - resVec[0].first)/stepSizes[k - 1];
  }	
	outRes = std::make_pair(outGrad, valid);
}

void TopOptObjFun::forwardDiffMPIDiff(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPIHandler* pMPIH, EvalFunction inef) const
{
	// default finite difference approximation
	std::vector<double> y;
	inTOR.getRealRep(y);
	std::vector<double> stepSizes(y.size());
	double sqrtEta = sqrt(getEps());
	std::cout << "Using MPI diff" << std::endl;
	for(std::size_t j = 0; j < y.size() ; ++j)
	{
		double stepSizeJ = sqrtEta*MAX(fabs(y[j]), inTOR.getDataMagnitude()) * (y[j] >= 0 ? 1. : -1.);
		if(y[j] == 1.) // Switch to backward difference, 1 is assumed max for all TORs
			stepSizeJ *= -1.;
		stepSizes[j] = stepSizeJ;
	}
	// MPI evaluate all TORs
	std::vector<std::pair<double, bool> > resVec;
	pMPIH->rootEvaluateDifference(y, stepSizes, resVec, inef);
	// Construct gradient from results
	std::vector<double> outGrad(inTOR.getDataSize());
	bool valid = false;
	if(!resVec.empty())
		valid = resVec[0].second;
	for(std::size_t k = 1; k < resVec.size(); ++k)
	{
		valid &= resVec[k].second;
		outGrad[k - 1] = (resVec[k].first - resVec[0].first)/stepSizes[k - 1];
	}
	outRes = std::make_pair(outGrad, valid);
}

void TopOptObjFun::centeredDiff(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPIHandler* pMPIH, EvalFunction inef) const
{
	std::vector<double> y;
	inTOR.getRealRep(y);
	std::vector<std::unique_ptr<TopOptRep> > torUPVec;
	torUPVec.reserve(2*y.size());
	std::vector<TopOptRep*> torVec(2*y.size());
	std::vector<double> x = y;
	std::vector<double> stepSizes(y.size());
	double sqrtEta = sqrt(getEps());
	for(std::size_t j = 0; j < y.size() ; ++j)
	{
		double stepSizeJ = sqrtEta*MAX(fabs(x[j]), inTOR.getDataMagnitude());
		// Forward half
		double tempJ = x[j];
		x[j] += stepSizeJ;
		double tempJp = x[j];
		torUPVec.push_back(inTOR.clone());
		torVec[2*j] = torUPVec.back().get();
		torVec[2*j]->setRealRep(x.begin(), x.end());
		// Backward half
		x[j] = tempJ - stepSizeJ;
		stepSizes[j] = tempJp - x[j];
		torUPVec.push_back(inTOR.clone());
		torVec[2*j + 1] = torUPVec.back().get();
		torVec[2*j + 1]->setRealRep(x.begin(), x.end());
		x[j] = tempJ;
	}
  // MPI evaluate all TORs
  std::vector<std::pair<double, bool> > resVec;
  pMPIH->rootBatchEvaluate(torVec, resVec, inef);
  // Construct gradient from results
  std::vector<double> outGrad(inTOR.getDataSize());
  bool valid = false;
  if(!resVec.empty())
    valid = resVec[0].second;
  for(std::size_t k = 0; k < stepSizes.size(); ++k)
  {
    valid &= resVec[k].second;
    outGrad[k] = (resVec[2*k].first - resVec[2*k + 1].first)/stepSizes[k];
  }
  outRes = std::make_pair(outGrad, valid);
}

void TopOptObjFun::gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	std::cout << "In serial gc" << std::endl;
	// Default finite difference approximation
	std::vector<double> y;
	inTOR.getRealRep(y);
	std::vector<double> x = y;
	double sqrtEta = sqrt(getEps());
	std::unique_ptr<TopOptRep> locTOR = inTOR.clone();
	std::pair<std::vector<double>, bool> cc;
	this->c(*locTOR, cc);
	std::vector<double> outGrad(y.size());
	bool valid = cc.second & (!cc.first.empty());
	for(std::size_t j = 0; j < y.size() && valid; ++j)
	{
		double stepSizeJ = sqrtEta*MAX(fabs(x[j]), inTOR.getDataMagnitude()) * (x[j] > 0 ? 1. : -1.);
		double tempJ = x[j];
		x[j] += stepSizeJ;
		stepSizeJ = x[j] - tempJ;
		locTOR->setRealRep(x.begin(), x.end());
		std::pair<std::vector<double>, bool> cj;
		this->c(*locTOR, cj);
		if(!cj.first.empty())
			outGrad[j] = (cj.first[0] - cc.first[0])/stepSizeJ;
		valid &= cj.second & (!cj.first.empty());
		x[j] = tempJ;
	}
	outRes = std::make_pair(outGrad, valid);
}

void TopOptObjFun::gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	computeGradient(efC, inTOR, outRes, communicator);
}

void TopOptObjFun::setEps()
{
	mEps = 1.;
	if (mEps == 1)
	{
		do
			mEps /= 2.;
		while (1. + mEps > 1.);
		mEps *= 2.;
	}
}

void TopOptObjFun::computeGradient(EvalFunction inEF, const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, 
	MPI::Comm& communicator) const
{
	outRes.first.clear();
	MPIHandler locMPIH(communicator, communicator.Get_size(), this, inTOR.clone());
	// Hold slave nodes
	if(locMPIH.amIRoot())
		forwardDiffMPIDiff(inTOR, outRes, &locMPIH, inEF);
	else
		locMPIH.slaveWaitAndProcessFlag();
}

void TopOptObjFun::fAndG(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& fRes,
												std::pair<std::vector<double>, bool>& gRes) const
{
	this->f(inTOR, fRes);
	this->g(inTOR, gRes);
}

void TopOptObjFun::fAndG(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& fRes,
													std::pair<std::vector<double>, bool>& gRes, MPI::Comm& communicator) const
{
	this->f(inTOR, fRes, communicator);
	this->g(inTOR, gRes, communicator);
}

}

