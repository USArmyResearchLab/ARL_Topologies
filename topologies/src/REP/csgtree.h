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

#ifndef CSGTREEREP_H
#define CSGTREEREP_H

#include "topoptrep.h"
#include "inputloaderrep.h"
#include "csgterminalnode.h"
#include <memory>

namespace Topologies{
class GenericMaterial;
class CSGNode;
class TOMesh;

//! A constructive solid geometry (CSG) topology optimization representation
/*! This implementation of TopOptRep is experimental and should be used with caution.
 *  It paramaterizes topology with a surface representation: A tree is used that defines
 *  a set of Boolean operations on solids defined as either the convex hull of a set of points
 *  or the alpha shape.  The optimization parameters are then the point coordinates.
 */
class CSGTreeRep : public TopOptRep
{
public:
	CSGTreeRep(const InputLoader::TORCSGTree& inputParams);
	CSGTreeRep(const std::vector<std::vector<int> >& discreteParams, const std::vector<std::vector<double> >& realParams);
	//! Constructs a 2D CSGTreeRep from a set of line segments
	/*! This will compute a convex partitioning of the polygons defined by the line segments and then create the 
	 *  tree from those polygons.
	 */
	CSGTreeRep(const InputLoader::TORCSGTree& inputParams, const std::vector<Mesh_Segment_2>& segVec);
	virtual ~CSGTreeRep();
	CSGTreeRep(const CSGTreeRep& copy);
	CSGTreeRep(CSGTreeRep && copy);
	CSGTreeRep& operator=(CSGTreeRep copy);
	void swap(CSGTreeRep& arg2);
	virtual std::unique_ptr<TopOptRep> clone() const;

	//! @name Decode functions
	//@{
	virtual void get2DSegments(std::vector<Mesh_Segment_2>& segVec) const;
	virtual std::unique_ptr<TOMesh> get2DMesh() const;
	virtual std::unique_ptr<TOMesh> get2DMesh(const GeometryTranslation::MesherData& meshParams) const;
	virtual void getBoundary(std::vector<Mesh_Segment_2>& boundaryVec) const;
	virtual std::unique_ptr<TOMesh> getOutputMesh() const;
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
	virtual void refine();
	virtual void prune();
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
	virtual void getDataSize(std::vector<std::size_t>& sizes) const;
	virtual double getDataMagnitude() const {return CSGTerminalNode::getBoundMag();}
	virtual unsigned getDimension() const;
	virtual double computeVolumeFraction() const;
	virtual std::vector<double> computeGradVolumeFraction() const {return std::vector<double>();}
	virtual std::vector<std::map<std::size_t, double>> diffRep() const {return std::vector<std::map<std::size_t, double>>();}
	virtual void getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
		std::vector<std::vector<double> >& realParams) const;
	virtual std::string getName() const {return getClassName();}
	//! Returns TopOptRep implementation's class name: csgtree
	static std::string getClassName() {return "csgtree";}
	//@}
	//! @name Data modification functions
	//@{
	virtual void filterData(std::vector<double>& valVec, double radius) const;
	virtual void filterData(double radius);
	virtual void boundsCheck(std::vector<double>& realVec) const;
	//@}
private:
	Nef_polyhedron_2 getNefPoly() const;
	Nef_polyhedron_2 getNefPolyWithBoundaries() const;
	void setMeshOptVals() const;
	void finishSetup();
	void finishSetup(const std::vector<Mesh_Segment_2>& segVec);
	GeometryTranslation::MesherData getDefaultMeshParams() const;
	void scaleParameters(std::vector<double>& newvals) const;
	void unscaleParameters(std::vector<double>& newvals) const;
	
	std::vector<double> affineVec;
	bool useAffine, shapesAreHoles;
	double width, height; // Region physical size
	unsigned nx, ny, nShapesX, nShapesY, nPointsPerShape;
	double minDensity;
	std::unique_ptr<CSGNode> rootNode;
	mutable CDT_2 mesh;
	static const double maxScale, maxMove;
private:
	// Hidden virtual functions
	virtual std::unique_ptr<TOMesh> get3DSurfaceMesh() const; // Define type later
	virtual std::unique_ptr<TOMesh> get3DVolumeMesh() const;
};

inline
	unsigned CSGTreeRep::getDimension() const
{
	return 2;
}
}
#endif

