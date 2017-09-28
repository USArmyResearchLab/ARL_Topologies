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

#ifndef TESTREP_H
#define TESTREP_H

#include "topoptrep.h"
#include "tomesh.h"
#include "helper.h"

namespace Topologies{
//! A mock TopOptRep to test the optimizers.  Uses 2 variables.
class TestRep : public TopOptRep
{
public:
	TestRep() : TopOptRep(tortPixel), optValx(0.), optValy(0.) {}
	TestRep(const std::vector<std::vector<int> >& discreteParams, const std::vector<std::vector<double> >& realParams) : 
		TopOptRep(tortPixel) {}
	virtual ~TestRep() {}
	virtual std::unique_ptr<TopOptRep> clone() const {return std::unique_ptr<TopOptRep>(new TestRep(*this));}

	//! @name Decode functions
	//@{
	virtual void get2DSegments(std::vector<Mesh_Segment_2>& segVec) const {}
	virtual std::unique_ptr<TOMesh> get2DMesh() const {return nullptr;}
	virtual std::unique_ptr<TOMesh> get2DMesh(const GeometryTranslation::MesherData& meshParams) const {return nullptr;}
	virtual void getBoundary(std::vector<Mesh_Segment_2>& boundaryVec) const {}
	virtual std::unique_ptr<TOMesh> getOutputMesh() const {return nullptr;}
	//@}
	//! @name Initialization functions
	//@{
	virtual void initialize() {optValx = optValy = 0.;}
	virtual void initialize(double val, std::pair<double, double> randRange) 
	{	optValx = val + HelperNS::RandomGen::instance().randRealInRange(randRange.first, randRange.second); 
		optValy = val + HelperNS::RandomGen::instance().randRealInRange(randRange.first, randRange.second);}
	virtual void initialize(double val) {optValx = optValy = val;}
	virtual void randomize() 
	{ optValx = HelperNS::RandomGen::instance().randRealInRange<double>();
		optValy = HelperNS::RandomGen::instance().randRealInRange<double>();}
	//@}
	//! @name Structural modification functions
	//@{
	virtual void refine() {}
	virtual void prune() {}
	//@}
	//! @name Data access
	//@{
	virtual void setRealRep(const std::vector<double>& newvals) 
	{assert(newvals.size() == 2); optValx = newvals[0]; optValy = newvals[1];}
	virtual void setDiscreteRep(const std::vector<int>& newvals) 
	{assert(newvals.size() == 2); optValx = (double)newvals[0]; optValy = (double)newvals[1];}
	virtual void setMPIRep(const std::vector<std::vector<int> >& discreteVars, const std::vector<std::vector<double> >& realVars)
	{ assert(!realVars.empty()); assert(realVars[0].size() == 2); optValx = realVars[0][0]; optValy = realVars[0][1];}
	virtual void getRealRep(std::vector<double>& realVec) const 
	{realVec = {optValx, optValy};}
	virtual void getDiscreteRep(std::vector<int>& discVec) const 
	{discVec = {(int)optValx, (int)optValy};}
	virtual void getMPIRep(std::vector<std::vector<int> >& discreteVars, std::vector<std::vector<double> >& realVars) const
	{realVars.resize(1); realVars[0] = {optValx, optValy};}
	virtual std::size_t getDataSize() const {return 2;}
	virtual void getDataSize(std::vector<std::size_t>& sizes) const {sizes = {2, 1};}
	virtual unsigned getDimension() const {return 2;}
	virtual double computeVolumeFraction() const {return optValx;}
	virtual std::vector<double> computeGradVolumeFraction() const {return std::vector<double>(2, 1.);}
	virtual std::vector<std::map<std::size_t, double>> diffRep() const {return std::vector<std::map<std::size_t, double>>();}
	virtual double getDataMagnitude() const {return 1.;}
	virtual void getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
		std::vector<std::vector<double> >& realParams) const {}
	virtual std::string getName() const {return getClassName();}
	//! Returns TopOptRep implementation's class name: testrep
	static std::string getClassName() {return "testrep";}
	//@}
	//! @name Data modification functions
	//@{
	virtual void filterData(std::vector<double>& valVec, double radius) const {}
	virtual void filterData(double radius) {}
	virtual void boundsCheck(std::vector<double>& realVec) const {}
	//@}

protected:
	double optValx, optValy;

private:
	// Hidden virtual functions
	virtual std::unique_ptr<TOMesh> get3DSurfaceMesh() const {return nullptr;} // Define type later
	virtual std::unique_ptr<TOMesh> get3DVolumeMesh() const {return nullptr;}
};
}
#endif

