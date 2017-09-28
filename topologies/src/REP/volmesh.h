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

#ifndef VOLMESH_H
#define VOLMESH_H

#include "topoptrep.h"
#include "inputloaderrep.h"
#include "helper.h"
#include "filterbase.h"
#include <memory>
#include <functional>

namespace Topologies{
class GenericMaterial;

//! Abstract base class for topology representations that discretize an arbitrary region with volume elments
/*! The class uses two template parameters which are a penalization function and a 
 *  projection function.  The penalization function should penalize intermediate values such as
 *  a SIMP or mod-SIMP (mod-SIMP is the default).  The projection function uses both a density
 *  filter and the specified projection to aid optimization.  Default is nothing, but it can be
 *  a Heaviside projection.
 */
template <typename PenaltyFunc = HelperNS::powPenalMin, typename ProjectionFunc = HelperNS::defaultProjFunc>
class VolMesh : public TopOptRep
{
public:
	explicit VolMesh(TORType inTORT) : TopOptRep(inTORT) {}
	VolMesh(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, const std::vector<double>& inPenalParams,
    const std::vector<double>& inProjParams);
	VolMesh(TORType inTORT, const InputLoader::TORGenericMesh& inputParams, const std::vector<double>& inPenalParams, 
		const std::vector<double>& inProjParams);
	virtual ~VolMesh(){}
	VolMesh(const VolMesh& copy);
	VolMesh(VolMesh && copy) : TopOptRep(copy.myTORT) {swap(copy);}
	void swap(VolMesh& arg2);
	virtual std::unique_ptr<TopOptRep> clone() const = 0;

	//! @name Decode functions
	//@{
	virtual void get2DSegments(std::vector<Mesh_Segment_2>& segVec) const = 0;
	virtual std::unique_ptr<TOMesh> get2DMesh() const = 0;
	virtual std::unique_ptr<TOMesh> get2DMesh(const GeometryTranslation::MesherData& meshParams) const = 0;
	virtual std::unique_ptr<TOMesh> get3DSurfaceMesh() const = 0;
	virtual std::unique_ptr<TOMesh> get3DVolumeMesh() const = 0;
	virtual void getBoundary(std::vector<Mesh_Segment_2>& boundaryVec) const = 0;
	virtual std::unique_ptr<TOMesh> getOutputMesh() const = 0;
	//@}
	//! @name Initialization functions
	//@{
	virtual void initialize();
	virtual void initialize(double val, std::pair<double, double> randRange);
	virtual void initialize(double val);
	virtual void randomize();
	//@}
	//! @name Structural modification functions
	//@{
	virtual void refine() = 0;
	virtual void prune() = 0;
	//@}
	//! @name Data access
	//@{
	virtual void setRealRep(const std::vector<double>& newvals);
	virtual void setDiscreteRep(const std::vector<int>& newvals);
	virtual void setMPIRep(const std::vector<std::vector<int> >& discreteVars, const std::vector<std::vector<double> >& realVars);
	virtual void getRealRep(std::vector<double>& realVec) const;
	virtual void getDiscreteRep(std::vector<int>& discVec) const;
	virtual void getMPIRep(std::vector<std::vector<int> >& discreteVars, std::vector<std::vector<double> >& realVars) const;
	virtual std::size_t getDataSize() const;
	virtual void getDataSize(std::vector<std::size_t>& sizes) const = 0;
	virtual unsigned getDimension() const = 0;
	virtual double computeVolumeFraction() const;
	virtual std::vector<double> computeGradVolumeFraction() const;
	virtual std::vector<std::map<std::size_t, double>> diffRep() const;
	virtual double getDataMagnitude() const {return 1.;}
	virtual void getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
		std::vector<std::vector<double> >& realParams) const = 0;
	virtual std::string getName() const {return getClassName();}
	//! Returns TopOptRep implementation's class name: volmesh
	static std::string getClassName() {return "volmesh";}
	//@}
	//! @name Data modification functions
	//@{
	virtual void filterData(std::vector<double>& valVec, double radius) const;
	virtual void filterData(double radius);
	virtual void boundsCheck(std::vector<double>& realVec) const;
	//@}

protected:
	void initFixedVals();
	void setFixedVals();
	void setFixedVals(std::vector<double>& inVec) const;
	std::vector<double> getFixedValElemMask() const;
	void boundsCheck();
	void setMeshOptVals(const std::vector<double>& optVals, TOMesh* pTOM) const;
	std::vector<double> getElemDensities() const;
	void computeDiffFilt();
	std::vector<double> getNodalAvgElemDensities() const;
	std::vector<std::map<std::size_t,double>> getDiffNodalAvgElemDensities() const;
	std::vector<std::map<std::size_t,double>> getIdentityFilt() const;
	std::vector<double> getDiffProjElemDensities() const;
	std::vector<double> getDiffProjElemDensities(const std::vector<double>& vals) const;
	std::vector<double> getProjElemDensities() const;
	std::vector<double> getProjElemDensities(const std::vector<double>& vals) const;
	std::vector<double> filterDensitiesWithDiffFilt() const;
	std::vector<Point_2_base> getElemCentroids2D() const;
	std::vector<Point_3_base> getElemCentroids3D() const;
protected:
	double threshold, filtRad; // threshold for generating discrete representation
	GeometryTranslation::MesherData meshParams;
	std::string fileName;
	VolMeshTORSpecification vmTORSpec;
	std::vector<double> penalParams, projParams;
	std::vector<std::pair<unsigned, double>> fixedBlockVec;

	std::vector<double> pixelArray, fixedVals;
	std::vector<std::map<std::size_t,double>> diffFilt;
	std::vector<std::size_t> fixedElemIDVec;
	mutable std::unique_ptr<TOMesh> upMesh;
	std::unique_ptr<FilterBase> upFilt, upDataFilt;
};

template <typename PenaltyFunc, typename ProjectionFunc>
VolMesh<PenaltyFunc, ProjectionFunc>::VolMesh(const VolMesh<PenaltyFunc, ProjectionFunc>& copy) :
	TopOptRep(copy),
	threshold(copy.threshold),
	filtRad(copy.filtRad),
	meshParams(copy.meshParams),
	fileName(copy.fileName),
	pixelArray(copy.pixelArray),
	diffFilt(copy.diffFilt),
	vmTORSpec(copy.vmTORSpec),
	penalParams(copy.penalParams),
	projParams(copy.projParams),
	fixedBlockVec(copy.fixedBlockVec),
	fixedElemIDVec(copy.fixedElemIDVec),
	upMesh(copy.upMesh->clone()),
	upFilt(copy.upFilt->clone()),
	upDataFilt(copy.upDataFilt ? copy.upDataFilt->clone() : nullptr)
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
void VolMesh<PenaltyFunc, ProjectionFunc>::swap(VolMesh<PenaltyFunc, ProjectionFunc>& arg2)
{
	TopOptRep::swap(arg2);
	std::swap(threshold, arg2.threshold);
	std::swap(filtRad, arg2.filtRad);
	std::swap(meshParams, arg2.meshParams);
	fileName.swap(arg2.fileName);
	pixelArray.swap(arg2.pixelArray);
	std::swap(diffFilt, arg2.diffFilt);
	std::swap(vmTORSpec, arg2.vmTORSpec);
	penalParams.swap(arg2.penalParams);
	projParams.swap(arg2.projParams);
	fixedBlockVec.swap(arg2.fixedBlockVec);
	fixedElemIDVec.swap(arg2.fixedElemIDVec);
	upMesh.swap(arg2.upMesh);
	upFilt.swap(arg2.upFilt);
	upDataFilt.swap(arg2.upDataFilt);
}
}//namespace
#endif

