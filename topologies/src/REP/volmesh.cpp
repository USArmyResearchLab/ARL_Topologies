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
using namespace HelperNS;

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh<PenaltyFunc, ProjectionFunc>::VolMesh(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, 
	const std::vector<double>& inPenalParams, const std::vector<double>& inProjParams) :
	TopOptRep(inTORT),
	threshold(inputParams.getThreshold()),
	filtRad(inputParams.getFiltRad()),
	vmTORSpec(inputParams.getVMTORS()),
	penalParams(inPenalParams),
	projParams(inProjParams),
	useInterpolatoryFilt(false)
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
	projParams(inProjParams),
	useInterpolatoryFilt(false)
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::initialize()
{
	std::size_t nelems = TOR::realOptVals.size();
	TOR::realOptVals = std::vector<double>(nelems, threshold);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::initialize(double val)
{
	std::size_t nelems = TOR::realOptVals.size();
	TOR::realOptVals = std::vector<double>(nelems, val);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::initialize(double val, std::pair<double, double> randRange)
{
	initialize(val);
	HelperNS::RGWrapper rgw(randRange);
	std::vector<double> noiseVec(TOR::realOptVals.size());
	std::generate(noiseVec.begin(), noiseVec.end(), rgw);
	std::transform(TOR::realOptVals.begin(), TOR::realOptVals.end(), noiseVec.begin(), TOR::realOptVals.begin(), std::plus<double>());
	boundsCheck();
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::randomize()
{
	HelperNS::RGWrapper rgw;
	std::generate(TOR::realOptVals.begin(), TOR::realOptVals.end(), rgw);
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
	setFixedVals(TOR::realOptVals);
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
	else if(vmTORSpec.torUnknownLocation == ulNode)
		return getNodalAvgElemDensities();
	return TOR::realOptVals;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::filterDensitiesWithDiffFilt() const
{
	return diffFilt.transposeTimes(TOR::realOptVals, upMesh->getNumElements());
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
void VolMesh<PenaltyFunc, ProjectionFunc>::computeDiffFilt(FilterBase const& filt)
{
	if(filtRad > 0.)
	{
		if(getDimension() == 2)
			diffFilt = filt.diffFilter(getElemCentroids2D(), filtRad);
		else
			diffFilt = filt.diffFilter(getElemCentroids3D(), filtRad);
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
			avg += TOR::realOptVals[*it];
		nodalAvgDens[k] = avg/(double)curElem.size();
	}
	return nodalAvgDens;
}

template <typename PenaltyFunc, typename ProjectionFunc>
HelperNS::SparseMatrix VolMesh<PenaltyFunc, ProjectionFunc>::getDiffNodalAvgElemDensities() const
{
	HelperNS::SparseMatrix res(upMesh->getNumNodes());
	for(std::size_t k = 0; k < upMesh->getNumElements(); ++k)
	{
		std::vector<std::size_t> curElem = upMesh->getElementConnectivity(k);
		for(auto it = curElem.begin(); it != curElem.end(); ++it)
			res.addEntry(*it, k, 1./static_cast<double>(curElem.size()));
	}
	return res;
}

template <typename PenaltyFunc, typename ProjectionFunc>
HelperNS::SparseMatrix VolMesh<PenaltyFunc, ProjectionFunc>::getIdentityFilt() const
{
	std::size_t fsize = vmTORSpec.torUnknownLocation == ulNode ? upMesh->getNumNodes() : upMesh->getNumElements();
	HelperNS::SparseMatrix res(fsize);
	for(std::size_t k = 0; k < fsize; ++k)
		res.addEntry(k, k, 1.);
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
void VolMesh<PenaltyFunc, ProjectionFunc>::updateRealRep()
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::setDiscreteRep(const std::vector<int>& newvals)
{
	TOR::realOptVals.resize(newvals.size());
	std::transform(newvals.begin(), newvals.end(), TOR::realOptVals.begin(), HelperNS::int2doub);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::setMPIRep(const std::vector<std::vector<int>>& discreteVars, 
	const std::vector<std::vector<double>>& realVars)
{
	if(!realVars.empty())
		setRealRep(realVars[0].begin(), realVars[0].end());
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::getRealRep(std::vector<double>& realVec) const
{
	realVec = TOR::realOptVals;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::getDiscreteRep(std::vector<int>& discVec) const
{
	discVec.resize(TOR::realOptVals.size());
	HelperNS::r2d locR2D(threshold);
	std::transform(TOR::realOptVals.begin(), TOR::realOptVals.end(), discVec.begin(), locR2D);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::getMPIRep(std::vector<std::vector<int>>& discreteVars, 
	std::vector<std::vector<double>>& realVars) const
{
	discreteVars.clear();
	realVars = {TOR::realOptVals};
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::size_t VolMesh<PenaltyFunc, ProjectionFunc>::getDataSize() const
{
	return TOR::realOptVals.size();
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
void VolMesh<PenaltyFunc, ProjectionFunc>::setupQuantitiesForDiffRep(std::vector<double>& rhoe, std::vector<double>& hprime, 
	std::vector<double>& fvmask, bool usePenalization) const
{
	// Get the filtered design variables
  std::vector<double> mue = getElemDensities();
	// Get the projection function of the filtered vals & its derivative
	hprime = getDiffProjElemDensities(mue);
	// Get fixed value mask, vector has 0's where fixed values are, i.e. don't contribute to gradient
	fvmask = getFixedValElemMask();
	// Get derivative of penalty function, if used
	if(usePenalization)
	{
		rhoe = getProjElemDensities(mue);
		PenaltyFunc dP(penalParams, true);
		std::transform(rhoe.begin(), rhoe.end(), rhoe.begin(), dP);
	}
	else
		rhoe = std::vector<double>(hprime.size(), 1.);
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
	SparseMatrix const& dF = diffFilt;
	// Get other quantities
	std::vector<double> hprime, fvmask, rhoe;
	setupQuantitiesForDiffRep(rhoe, hprime, fvmask, false);
	// Gradient of volume fraction is then just the sum of each row
	std::vector<double> grad(dF.size(),0.);
	for(std::size_t k = 0; k < grad.size(); ++k)
	{
		SparseMatrix::SparseRow const& curRow = dF.row(k);
		for(auto const & colPair : curRow)
		{
			std::size_t kcol = SparseMatrix::index(colPair);
			grad[k] += elemAreas[kcol]*hprime[kcol]*fvmask[kcol]*SparseMatrix::value(colPair)/totArea;
		}
	}
	return grad;
}

template <typename PenaltyFunc, typename ProjectionFunc>
std::vector<double> VolMesh<PenaltyFunc, ProjectionFunc>::applyDiffRep(std::vector<double> const& elemGrad, 
	bool usePenalization) const
{
	// First get sparse matrix of partial derivatives of filter wrt the design variables
	SparseMatrix const& dF = diffFilt;
	// Get other quantities
	std::vector<double> hprime, fvmask, rhoe;
	setupQuantitiesForDiffRep(rhoe, hprime, fvmask, usePenalization);
	// Form full set of partial derivatives
	std::vector<double> fullGrad(dF.size(),0.);
	auto gradIt = fullGrad.begin();
	for(auto& curRow : dF)
	{
		for(auto& colPair : curRow)
		{
			std::size_t kcol = SparseMatrix::index(colPair);
			*gradIt += elemGrad[kcol]*fvmask[kcol]*hprime[kcol]*rhoe[kcol]*SparseMatrix::value(colPair);
		}
		++gradIt;
	}
	return fullGrad;
}

template <typename PenaltyFunc, typename ProjectionFunc>
HelperNS::SparseMatrix VolMesh<PenaltyFunc, ProjectionFunc>::diffRep(bool usePenalization) const
{
	// First get sparse matrix of partial derivatives of filter wrt the design variables
	SparseMatrix dF = diffFilt;
	// Get other quantities
	std::vector<double> hprime, fvmask, rhoe;
	setupQuantitiesForDiffRep(rhoe, hprime, fvmask, usePenalization);
	// Form full set of partial derivatives
	for(auto& curRow : dF)
	{
		for(auto& colPair : curRow)
		{
			std::size_t kcol = SparseMatrix::index(colPair);
			SparseMatrix::value(colPair) *= fvmask[kcol]*hprime[kcol]*rhoe[kcol];
		}
	}
	return dF;
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::filterData(std::vector<double>& valVec, double radius) const
{
	if(!upDataFilt) // Only construct if needed
		upDataFilt = constructFilter();
	valVec = (*upDataFilt)(valVec, radius);
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::filterData(double radius)
{
	filterData(TOR::realOptVals, radius);
}

template class VolMesh<HelperNS::powPenal, HelperNS::defaultProjFunc>;
template class VolMesh<HelperNS::powPenalMin, HelperNS::defaultProjFunc>;
template class VolMesh<HelperNS::powPenal, HelperNS::regularizedHeaviside>;
template class VolMesh<HelperNS::powPenalMin, HelperNS::regularizedHeaviside>;
template class VolMesh<HelperNS::powPenal, HelperNS::thresholdHeaviside>;
template class VolMesh<HelperNS::powPenalMin, HelperNS::thresholdHeaviside>;
} //namespace 

