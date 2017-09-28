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

#ifndef TOTESTOBJFUN2_H
#define TOTESTOBJFUN2_H

#include "topoptobjfun.h"
#include "topoptrep.h"

namespace Topologies{
//! A test class for testing the finite difference gradient of the base class TopOptObjFun
class TOTestObjFun2 : public TopOptObjFun
{
public:
	//! Constructor which takes the name of an input file.
	TOTestObjFun2(double inTarget, double inConstraint, bool inXconstraint = true) : 
		target(inTarget), constraint(inConstraint), xconstraint(inXconstraint) {}
	virtual ~TOTestObjFun2() {}

	virtual std::pair<double, bool> operator()(const TopOptRep& inTOR) const;
	virtual void f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
	virtual void c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
	virtual void g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPIHandler* pMPIH) const 
	{g(inTOR, outRes);}
  virtual void g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;

  virtual void gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const;
  virtual void gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPIHandler* pMPIH) const 
	{gc(inTOR, outRes);}
#ifdef USE_MPI
	virtual std::pair<double, bool> operator()(const TopOptRep& inTOR, MPI::Comm& communicator) const;
	virtual void f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
	virtual void c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const;
#endif
	virtual void printResult(const TopOptRep& inTOR, const std::string& fileName) const {}
	bool gradCheck(const TopOptRep& inTOR) const;

private:
	bool xconstraint;
	double target, constraint;
};

inline
std::pair<double, bool> TOTestObjFun2::operator()(const TopOptRep& inTOR) const
{
	std::pair<double, bool> retPair;
	std::vector<double> val;
	inTOR.getRealRep(val);
	retPair.first = (val[0] - target)*(val[0] - target) + (val[1] - target)*(val[1] - target);
	retPair.second = true;
	return retPair;
}

inline
void TOTestObjFun2::f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	std::pair<double, bool> res = (*this)(inTOR);
	outRes.first = {res.first};
	outRes.second = res.second;
}

inline
void TOTestObjFun2::c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	std::vector<double> val;
	inTOR.getRealRep(val);
	if(xconstraint)
		outRes.first = {val[0] - constraint};
	else
		outRes.first = {val[0] + val[1] - constraint};
	outRes.second = true;
}

inline
void TOTestObjFun2::g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	std::vector<double> val;
	inTOR.getRealRep(val);
	outRes.first = {2.*(val[0] - target), 2.*(val[1] - target)};
	outRes.second = true;
}

inline
void TOTestObjFun2::gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	if(xconstraint)
		outRes.first = {1., 0.};
	else
		outRes.first = {1., 1.};
	outRes.second = true;
}

inline
bool TOTestObjFun2::gradCheck(const TopOptRep& inTOR) const
{
	std::pair<std::vector<double>, bool> outRes1, outRes2;
	g(inTOR, outRes1);
	TopOptObjFun::g(inTOR, outRes2);
	bool res = outRes1.second && outRes2.second;
	for(std::size_t k = 0; k < outRes1.first.size(); ++k)
		res &= TOL_EQ(outRes1.first[k], outRes2.first[k], 1e-7);
	// Constraint
	gc(inTOR, outRes1);
	TopOptObjFun::gc(inTOR, outRes2);
	res &= outRes1.second && outRes2.second;
	for(std::size_t k = 0; k < outRes1.first.size(); ++k)
		res &= TOL_EQ(outRes1.first[k], outRes2.first[k], 1e-7);
	return res;
}

#ifdef USE_MPI
inline
std::pair<double, bool> TOTestObjFun2::operator()(const TopOptRep& inTOR, MPI::Comm& communicator) const
{
	return (*this)(inTOR);
}

inline
void TOTestObjFun2::f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	f(inTOR, outRes);
}

inline
void TOTestObjFun2::c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	c(inTOR, outRes);
}
#endif
}

extern "C" Topologies::TopOptObjFun* create(const std::string& inpFileName)
{
	return new Topologies::TOTestObjFun2(0., 0.);
}

extern "C" void destroy(Topologies::TopOptObjFun* p)
{
	delete p;
}

#endif

