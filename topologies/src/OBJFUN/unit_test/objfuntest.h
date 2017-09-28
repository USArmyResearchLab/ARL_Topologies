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

#ifndef TOTESTOBJFUN_H
#define TOTESTOBJFUN_H

#include "topoptobjfun.h"
#include "topoptrep.h"

namespace Topologies{
//! A test class for testing the finite difference gradient of the base class TopOptObjFun
class TOTestObjFun : public TopOptObjFun
{
public:
	//! Constructor which takes the name of an input file.
	TOTestObjFun() {}
	virtual ~TOTestObjFun() {}

	virtual std::pair<double, bool> operator()(const TopOptRep& inTOR) const;
	virtual void f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
	virtual void c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
#ifdef USE_MPI
	virtual std::pair<double, bool> operator()(const TopOptRep& inTOR, MPI::Comm& communicator) const;
	virtual void f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
	virtual void c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
#endif
	virtual void printResult(const TopOptRep& inTOR, const std::string& fileName) const;
};

inline
std::pair<double, bool> TOTestObjFun::operator()(const TopOptRep& inTOR) const
{
	std::pair<double, bool> retPair;
	retPair.first = inTOR.computeVolumeFraction();
	retPair.second = true;
	return retPair;
}

inline
void TOTestObjFun::f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	outRes.first.resize(1);
	std::pair<double, bool> res = (*this)(inTOR);
	outRes.first[0] = res.first;
	outRes.second = res.second;
}

inline
void TOTestObjFun::c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	f(inTOR, outRes);
}

#ifdef USE_MPI
inline
std::pair<double, bool> TOTestObjFun::operator()(const TopOptRep& inTOR, MPI::Comm& communicator) const
{
	return (*this)(inTOR);
}

inline
void TOTestObjFun::f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	f(inTOR, outRes);
}

inline
void TOTestObjFun::c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	c(inTOR, outRes);
}
#endif

inline
void TOTestObjFun::printResult(const TopOptRep& inTOR, const std::string& fileName) const 
{
	std::cout << "In printResult, got name: " << fileName << std::endl;
}
}//namespace

extern "C" Topologies::TopOptObjFun* create(const std::string& inpFileName)
{
	return new Topologies::TOTestObjFun();
}

extern "C" void destroy(Topologies::TopOptObjFun* p)
{
	delete p;
}

#endif
