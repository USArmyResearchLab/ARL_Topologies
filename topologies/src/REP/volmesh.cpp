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

#include "volmesh.h"
#include "geometrytranslation.h"
#include "helper.h"
#include "tomesh.h"
#include <algorithm>
#include <unordered_set>
#include "tomeshprocessing.h"

namespace Topologies{
template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh<PenaltyFunc, ProjectionFunc>::VolMesh(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, 
	const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
	TopOptRep(inTORT),
	threshold(inputParams.getThreshold()),
	filtRad(inputParams.getFiltRad()),
	vmTORSpec(inputParams.getVMTORS()),
	penalParams(inPenalParams),
	projParams(inProjParams)
{
}


template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh<PenaltyFunc, ProjectionFunc>::VolMesh(TORType inTORT, const InputLoader::TORGenericMesh& inputParams, 
	const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
	TopOptRep(inTORT),
	threshold(inputParams.getThreshold()),
	filtRad(inputParams.getFiltRad()),
	meshParams(inputParams.getMeshParams()),
	fileName(inputParams.getFileName()),
	fixedBlockVec(inputParams.getFixedBlockVec()),
	vmTORSpec(inputParams.getVMTORS()),
	penalParams(inPenalParams),
	projParams(inProjParams)
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::initialize()
{
	std::size_t nelems = pixelArray.size();
	pixelArray = std::vector<double>(nelems, threshold);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::initialize(double val)
{
	std::size_t nelems = pixelArray.size();
	pixelArray = std::vector<double>(nelems, val);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::initialize(double val, std::pair<double, double> randRange)
{
	initialize(val);
	HelperNS::RGWrapper rgw(randRange);
	std::vector<double> noiseVec(pixelArray.size());
	std::generate(noiseVec.begin(), noiseVec.end(), rgw);
	std::transform(pixelArray.begin(), pixelArray.end(), noiseVec.begin(), pixelArray.begin(), std::plus<double>());
	boundsCheck();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::randomize()
{
	HelperNS::RGWrapper rgw;
	std::generate(pixelArray.begin(), pixelArray.end(), rgw);
	boundsCheck();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::initFixedVals()
{
	if(fixedBlockVec.empty())
		return;
	// Set all fixed element densities and generate a vector of unknown ids
	fixedElemIDVec.clear();
	fixedElemIDVec.reserve(upMesh->getNumElements());
	for(std::size_t ke = 0; ke < upMesh->getNumElements(); ++ke)
	{
		bool free = true;
		unsigned matid = upMesh->getMatID(ke);
		for(std::size_t kfb = 0; kfb < fixedBlockVec.size() && free; ++kfb)
		{
			if(matid == fixedBlockVec[kfb].first)
			{
				free = false;
				upMesh->setOptVal(ke, fixedBlockVec[kfb].second);
				fixedElemIDVec.push_back(ke);
			}
		}
	}
	fixedElemIDVec.shrink_to_fit();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::setFixedVals()
{
	setFixedVals(pixelArray);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::setFixedVals(std::vector<double>& inVec) const
{
	if(fixedBlockVec.empty())
		return;
	for(std::size_t k = 0; k < fixedElemIDVec.size(); ++k)
	{
		if(fixedElemIDVec[k] >= inVec.size()) // Check for validity
			continue;
		// Get block id
		unsigned mid = upMesh->getMatID(fixedElemIDVec[k]);
		for(std::size_t kfb = 0; kfb < fixedBlockVec.size(); ++kfb)
		{
			if(fixedBlockVec[kfb].first == mid)
				inVec[fixedElemIDVec[k]] = fixedBlockVec[kfb].second; // Set value to fixed block value
		}
	}
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::getElemDensities() const
{
	if(filtRad > 0.)
		return filterDensitiesWithDiffFilt();
/*	{
		if(getDimension() == 2)
			return filterDensitiesWithDiffFilt();
			return (*upFilt)(pixelArray, getElemCentroids2D(), filtRad);
		else
			return (*upFilt)(pixelArray, getElemCentroids3D(), filtRad);
	}*/
	else if(vmTORSpec.torUnknownLocation == ulNode)
		return getNodalAvgElemDensities();
	return pixelArray;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::filterDensitiesWithDiffFilt() const
{
	std::vector<double> res(upMesh->getNumElements(), 0.);
	for(std::size_t krow = 0; krow < diffFilt.size(); ++krow)
	{
		const std::map<std::size_t,double>& curRow = diffFilt[krow];
		for(auto colIt = curRow.begin(); colIt != curRow.end(); ++colIt)
			res[colIt->first] += pixelArray[krow]*(colIt->second);
	}
	return res;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<Point_2_base> VolMesh<PenaltyFunc, ProjectionFunc>::getElemCentroids2D() const
{
	std::vector<Point_2_base> ptVec(upMesh->getNumElements());
	for(std::size_t k = 0; k < upMesh->getNumElements(); ++k)
		ptVec[k] = TOMeshProcessing::getElementCentroid2D(k, upMesh.get());
	return ptVec;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<Point_3_base> VolMesh<PenaltyFunc, ProjectionFunc>::getElemCentroids3D() const
{
	std::vector<Point_3_base> ptVec(upMesh->getNumElements());
	for(std::size_t k = 0; k < upMesh->getNumElements(); ++k)
		ptVec[k] = TOMeshProcessing::getElementCentroid3D(k, upMesh.get());
	return ptVec;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::computeDiffFilt()
{
	if(filtRad > 0.)
	{
		if(getDimension() == 2)
			diffFilt = upFilt->diffFilter(getElemCentroids2D(), filtRad);
		else
			diffFilt = upFilt->diffFilter(getElemCentroids3D(), filtRad);
	}
	else if(vmTORSpec.torUnknownLocation == ulNode)
		diffFilt = getDiffNodalAvgElemDensities();
	else
		diffFilt = getIdentityFilt();
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::getProjElemDensities() const
{
	return getProjElemDensities(getElemDensities());
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::getProjElemDensities(const std::vector<double>& vals) const
{
	ProjectionFunc projf(projParams);
	std::vector<double> res(vals.size());
	std::transform(vals.begin(), vals.end(), res.begin(), projf);
	return res;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::getDiffProjElemDensities() const
{
	return getDiffProjElemDensities(getElemDensities());
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::getDiffProjElemDensities(const std::vector<double>& vals) const
{
	ProjectionFunc dH(projParams, true);
	std::vector<double> res(vals.size());
	std::transform(vals.begin(), vals.end(), res.begin(), dH);
	return res;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::getNodalAvgElemDensities() const
{
	std::vector<double> nodalAvgDens(upMesh->getNumElements());
	for(std::size_t k = 0; k < nodalAvgDens.size(); ++k)
	{
		std::vector<std::size_t> curElem = upMesh->getElementConnectivity(k);
		double avg = 0.;
		for(auto it = curElem.begin(); it != curElem.end(); ++it)
			avg += pixelArray[*it];
		nodalAvgDens[k] = avg/(double)curElem.size();
	}
	return nodalAvgDens;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<std::map<std::size_t,double>> VolMesh<PenaltyFunc, ProjectionFunc>::getDiffNodalAvgElemDensities() const
{
	std::vector<std::map<std::size_t,double>> res(upMesh->getNumNodes());
	for(std::size_t k = 0; k < upMesh->getNumElements(); ++k)
	{
		std::vector<std::size_t> curElem = upMesh->getElementConnectivity(k);
		for(auto it = curElem.begin(); it != curElem.end(); ++it)
		{
			std::map<std::size_t, double>& curRow = res[*it];
			curRow[k] = 1./(double)curElem.size();
		}
	}
	return res;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<std::map<std::size_t,double>> VolMesh<PenaltyFunc, ProjectionFunc>::getIdentityFilt() const
{
	std::size_t fsize = vmTORSpec.torUnknownLocation == ulNode ? upMesh->getNumNodes() : upMesh->getNumElements();
	std::vector<std::map<std::size_t,double>> res(fsize);
	for(std::size_t k = 0; k < fsize; ++k)
		res[k][k] = 1.;
	return res;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::setMeshOptVals(const std::vector<double>& optVals, TOMesh* pTOM) const
{
	if(optVals.size() == pTOM->getNumElements())
		pTOM->setOptVals(optVals);
}

// Data access
template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::setRealRep(const std::vector<double>& newvals)
{
	assert(newvals.size() == pixelArray.size());
	pixelArray = newvals;
	boundsCheck();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::setDiscreteRep(const std::vector<int>& newvals)
{
	pixelArray.resize(newvals.size());
	std::transform(newvals.begin(), newvals.end(), pixelArray.begin(), HelperNS::int2doub);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::setMPIRep(const std::vector<std::vector<int> >& discreteVars, 
	const std::vector<std::vector<double> >& realVars)
{
	if(!realVars.empty())
		setRealRep(realVars[0]);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::getRealRep(std::vector<double>& realVec) const
{
	realVec = pixelArray;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::getDiscreteRep(std::vector<int>& discVec) const
{
	discVec.resize(pixelArray.size());
	HelperNS::r2d locR2D(threshold);
	std::transform(pixelArray.begin(), pixelArray.end(), discVec.begin(), locR2D);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::getMPIRep(std::vector<std::vector<int> >& discreteVars, 
	std::vector<std::vector<double> >& realVars) const
{
	discreteVars.clear();
	realVars = {pixelArray};
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::size_t VolMesh<PenaltyFunc, ProjectionFunc>::getDataSize() const
{
	return pixelArray.size();
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::getFixedValElemMask() const
{
	std::vector<double> fvmask(upMesh->getNumElements(), 1.);
	for(std::size_t k = 0; k < fixedElemIDVec.size(); ++k)
		fvmask[fixedElemIDVec[k]] = 0.;
	return fvmask;
}

template <typename PenaltyFunc, typename ProjectionFunc>
double VolMesh<PenaltyFunc, ProjectionFunc>::computeVolumeFraction() const
{
	std::vector<double> projDensities = getProjElemDensities();
	setFixedVals(projDensities);
	double propTot = 0., totArea = 0.;
	for(std::size_t k = 0; k < upMesh->getNumElements(); ++k)
	{
		Real area = TOMeshProcessing::computeElementVolume(k, upMesh.get());
		propTot += projDensities[k]*area;
		totArea += area;
	}
	return propTot/totArea;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::computeGradVolumeFraction() const
{
	// Compute area
	std::vector<double> elemAreas(upMesh->getNumElements());
	for(std::size_t k = 0; k < upMesh->getNumElements(); ++k)
  	elemAreas[k] = TOMeshProcessing::computeElementVolume(k, upMesh.get());
	double totArea = std::accumulate(elemAreas.begin(), elemAreas.end(), 0.);
	// First get sparse matrix of partial derivatives of filter wrt the design variables
	std::vector<std::map<std::size_t,double>> dF = diffFilt;
	// Next, get the filtered design variables
	std::vector<double> mue = getElemDensities();
	// Get the projection function of the filtered vals & its derivative
	std::vector<double> hprime = getDiffProjElemDensities(mue);
	// Get fixed value mask, vector has 0's where fixed values are, i.e. don't contribute to gradient
	std::vector<double> fvmask = getFixedValElemMask();
	// Gradient of volume fraction is then just the sum of each row
	std::vector<double> grad(dF.size(),0.);
	for(std::size_t k = 0; k < grad.size(); ++k)
	{
		std::map<std::size_t,double>& curRow = dF[k];
		for(auto it = curRow.begin(); it != curRow.end(); ++it)
			grad[k] += elemAreas[it->first]*hprime[it->first]*fvmask[it->first]*(it->second)/totArea;
	}
	return grad;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<std::map<std::size_t, double>> VolMesh<PenaltyFunc, ProjectionFunc>::diffRep() const
{
	// First get sparse matrix of partial derivatives of filter wrt the design variables
	std::vector<std::map<std::size_t,double>> dF = diffFilt;
	// Next, get the filtered design variables
	std::vector<double> mue = getElemDensities();
	// Get the projection function of the filtered vals & its derivative
	std::vector<double> hprime = getDiffProjElemDensities(mue);
	std::vector<double> rhoe = getProjElemDensities(mue);
	// Get derivative of penalty function
	PenaltyFunc dP(penalParams, true);
	std::transform(rhoe.begin(), rhoe.end(), rhoe.begin(), dP);	
	// Get fixed value mask, vector has 0's where fixed values are, i.e. don't contribute to gradient
	std::vector<double> fvmask = getFixedValElemMask();
	// Form full set of partial derivatives
	for(auto it = dF.begin(); it != dF.end(); ++it)
	{
		std::map<std::size_t,double>& curRow = *it;
		for(auto mit = curRow.begin(); mit != curRow.end(); ++mit)
			mit->second = fvmask[mit->first]*(mit->second)*hprime[mit->first]*rhoe[mit->first];
	}
	return dF;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::boundsCheck(std::vector<double>& realVec) const
{
	std::replace_if(realVec.begin(), realVec.end(), HelperNS::greaterThan1, 1.);
	std::replace_if(realVec.begin(), realVec.end(), HelperNS::lessThan0, 0.);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::boundsCheck()
{
	std::replace_if(pixelArray.begin(), pixelArray.end(), HelperNS::greaterThan1, 1.);
	std::replace_if(pixelArray.begin(), pixelArray.end(), HelperNS::lessThan0, 0.);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::filterData(std::vector<double>& valVec, double radius) const
{
	if(upDataFilt)
		valVec = (*upDataFilt)(valVec, radius);
	else
		valVec = (*upFilt)(valVec, radius);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::filterData(double radius)
{
	filterData(pixelArray, radius);
}

template class VolMesh<HelperNS::powPenal, HelperNS::defaultProjFunc>;
template class VolMesh<HelperNS::powPenalMin, HelperNS::defaultProjFunc>;
template class VolMesh<HelperNS::powPenal, HelperNS::regularizedHeaviside>;
template class VolMesh<HelperNS::powPenalMin, HelperNS::regularizedHeaviside>;
template class VolMesh<HelperNS::powPenal, HelperNS::thresholdHeaviside>;
template class VolMesh<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>;
} //namespace 
