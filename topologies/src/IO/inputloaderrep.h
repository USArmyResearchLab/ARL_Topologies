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

#ifndef INPUTLOADERREP_H
#define INPUTLOADERREP_H

#include "inputloader.h"

namespace Topologies{
namespace InputLoader
{
	//! Class to handle loading of PixelRep and VoxelRep
	class TORGenericVolume : public InputParser
	{
	public:
		TORGenericVolume(const std::string& inRepName) : repName(inRepName) {}
		virtual ~TORGenericVolume() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns the representation name
		const std::string& getRepName() const {return repName;}
		//! Returns the threshold for post-processing
		double getThreshold() const {return threshold;}
		//! Returns the analysis penalty power
		double getPenalPower() const {return penalPower;}
		//! Returns the minimum element denstity
		double getMinDensity() const {return minDensity;}
		//! Returns the density filter radius
		double getFiltRad() const {return filtRad;}
		//! Returns the regularized Heaviside parameter
		double getBetaHeavi() const {return betaHeavi;}
		//! Returns the region physical dimensions
		double getDimensions(std::size_t k) const {assert(k < dimensions.size()); return dimensions[k];}
		//! Returns the discretization sizes
		unsigned getDiscSizes(std::size_t k) const {assert(k < discSizes.size()); return discSizes[k];}
		//! Returns the rank for LowRankPixel
		unsigned getRank() const {return rank;}
		//! Returns the MeshElementType
		MeshElementType getMET() const {return myMET;}
		//! Returns the VolMeshTORSpecification object 
		const VolMeshTORSpecification& getVMTORS() const {return myVMTORS;}
	private:
		void setMeshElementType(const pugi::xml_node& rootNode);

		const std::string repName;
		double threshold = 0.5;
		double penalPower = 1., minDensity = 0., filtRad = 0., betaHeavi = 1e-8;
		std::vector<double> dimensions;
		std::vector<unsigned> discSizes;
		unsigned rank = 10;
		MeshElementType myMET;
		VolMeshTORSpecification myVMTORS;
	};

	class TORCSGTree : public InputParser
	{
	public:
		TORCSGTree(const std::string& inRepName) : repName(inRepName) {}
		virtual ~TORCSGTree() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns the representation name
		std::string getRepName() const {return repName;}
		//! Returns the region physical size
		double getRegionDimensions(std::size_t k) const {assert(k < regionDimensions.size()); return regionDimensions[k];}
		//! Returns the mesh discretization sizes
		unsigned getMeshSizes(std::size_t k) const {assert(k < meshSizes.size()); return meshSizes[k];}
		//! Returns the number of convex shapes
		unsigned getShapeNums(std::size_t k) const {assert(k < shapeNums.size()); return shapeNums[k];}
		//! Returns the number of points per shape
		unsigned getNumPointsPerShape() const {return nPointsPerShape;}
		//! Returns the minimum element density
		double getMinDensity() const {return minDensity;}
		//! Returns whether or not to use affine transformations as the optimization parameters
		bool getUseAffine() const {return useAffine;}
		//! Returns whether or not the CSG shapes represent holes
		bool getShapesAreHoles() const {return shapesAreHoles;}
	private:
		std::string repName;
		std::vector<double> regionDimensions;
		std::vector<unsigned> meshSizes, shapeNums;
		unsigned nPointsPerShape;
		double minDensity = 0.;
		bool useAffine = false, shapesAreHoles = false;
	};

	class CGALMesh : public InputParser
	{
	public:
		CGALMesh() {}
		virtual ~CGALMesh() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns the MesherData
		const GeometryTranslation::MesherData& getMesherData() const {return inputParams;}
	private:
		GeometryTranslation::MesherData inputParams;
		unsigned dim;
	};

	class TORGenericMesh : public InputParser
	{
	public:
		TORGenericMesh(const std::string& inRepName) : repName(inRepName) {}
		virtual ~TORGenericMesh() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns the representation name
		std::string getRepName() const {return repName;}
		//! Returns a vector of polygons to mesh
		const std::vector<std::vector<Point_2_base> >& getPolyVec() const {return ptVecs;}
		//! Returns the threshold for post-processing
		double getThreshold() const {return threshold;}
		//! Returns the analysis penalty power
		double getPenalPower() const {return penalPower;}
		//! Returns the minimum density
		double getMinDensity() const {return minDensity;}
		//! Returns the density filter radius
		double getFiltRad() const {return filtRad;}
		//! Returns the regularized Heaviside parameter
		double getBetaHeavi() const {return betaHeavi;}
		//! Returns the CGALMesh parameters
		const GeometryTranslation::MesherData& getMeshParams() const {return meshParams.getMesherData();}
		//! Returns the mesh file name
		const std::string& getFileName() const {return fileName;}
		//! Returns the vector of fixed block definitions
		const std::vector<std::pair<unsigned, double>>& getFixedBlockVec() const {return fixedBlockVec;}
		//! Returns the VolMeshTORSpecification object
		const VolMeshTORSpecification& getVMTORS() const {return myVMTORS;}
	private:
		void setFileFormat(const pugi::xml_node& rootNode);
		void setPolygonList(const pugi::xml_node& rootNode);
		std::vector<Point_2_base> parsePointVec(const pugi::xml_node& rootNode);

		std::string repName;
		std::vector<std::vector<Point_2_base> > ptVecs;
		double threshold = 0.5, penalPower = 1., minDensity = 0., filtRad = 0., betaHeavi = 1e-8;;
		CGALMesh meshParams;
		std::string fileName;
		std::vector<std::pair<unsigned, double>> fixedBlockVec;
		VolMeshTORSpecification myVMTORS;
	};
}
}
#endif

