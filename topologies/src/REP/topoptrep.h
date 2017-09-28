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

#ifndef TOPOPTREP_H
#define TOPOPTREP_H

#include <memory>
#include <string>
#include "cgal_types.h"
#include "geometrytranslation.h"
#include "topologiesdefs.h"

namespace Topologies{
class TOMesh;

//! Abstract base class for a topology optimization representation
/*! The main purpose of this class is to parameterize a topology using some given method (pixel, CSG, alpha shape, etc.)
 *  so that those parameters can be manipulated by an optimization algorithm.  The parameters are used to map the topology
 *  to a mesh, which is used for analysis.
 *  Implementations need to be able to return a list of segments or tri mesh for 2D or facets or tet mesh for 3D.
 */
class TopOptRep
{
public:
	virtual ~TopOptRep(){}
	explicit TopOptRep(TORType inTORT) : myTORT(inTORT) {}
	TopOptRep(const TopOptRep& copyFrom) : myTORT(copyFrom.myTORT) {}
	TopOptRep& operator=(TopOptRep copy) {swap(copy); return *this;}
	TopOptRep(TopOptRep && rhs) {swap(rhs);}
	void swap(TopOptRep& rhs) {std::swap(myTORT, rhs.myTORT);}
	virtual std::unique_ptr<TopOptRep> clone() const = 0;

	//! @name Decode functions
	//@{
	//! Generates a set of 2D line segments representing the topology
	virtual void get2DSegments(std::vector<Mesh_Segment_2>& segVec) const = 0;
	//! Returns a 2D mesh representing the topology to be used for analysis
	virtual std::unique_ptr<TOMesh> get2DMesh() const = 0;
	//! Returns a 2D mesh representing the topology with given mesh parameters to be used for analysis
	virtual std::unique_ptr<TOMesh> get2DMesh(const GeometryTranslation::MesherData& meshParams) const = 0;
	//! Returns a 3D surface mesh representing the topology to be used for analysis
	virtual std::unique_ptr<TOMesh> get3DSurfaceMesh() const = 0;
	//! Returns a 3D volume mesh representing the topology to be used for analysis
	virtual std::unique_ptr<TOMesh> get3DVolumeMesh() const = 0;
	//! Returns the boundary of the topology
	virtual void getBoundary(std::vector<Mesh_Segment_2>& boundaryVec) const = 0;
	//! Returns a mesh for output purposes (rather than analysis)
	virtual std::unique_ptr<TOMesh> getOutputMesh() const = 0;
	//@}
	//! @name Initialization functions
	//@{
	//! Initialize the TopOptRep to a default value
	virtual void initialize() = 0;
	//! Initialize all optimization parameters to the given value
	virtual void initialize(double val) = 0;
	//! Initialize all optimization parameters to random values centered at `val` and within `randRange`
	virtual void initialize(double val, std::pair<double, double> randRange) = 0;
	//! Randomize optimization parameters
	virtual void randomize() = 0;
	//@}
	//! @name Structural modification functions
	//@{
	//! Refine the structure of the topology so that it is higher fidelity, potentially with more optimization parameters
	virtual void refine() = 0;
	//! Corsen the structure of the topology so that it is lower fidelity, potentially with fewer optimization parameters
	virtual void prune() = 0;
	//@}
	//! @name Data access
	//@{
	//! Replaces the continuous optimization parameters of a topology with `newvals`
	/*! This is the main way that optimization algorithms (TopOpt) interact with a TopOptRep.  
	 *  Also, these values should always be between 0 and 1.
	 */
	virtual void setRealRep(const std::vector<double>& newvals) = 0;
	//! Replaces the discrete optimization parameters of a topology with `newvals`
	virtual void setDiscreteRep(const std::vector<int>& newvals) = 0;
	//! Replaces all optimization parameters with the input arguments
	/*! The values contained in the arguments must fully specify the topology.  This function is used to
	 *  receive values over MPI communication for parallel computations.
	 */
	virtual void setMPIRep(const std::vector<std::vector<int> >& discreteVars, const std::vector<std::vector<double> >& realVars) = 0;
	//! Copy the continuous topology representation into `realVec`
	/*! This is the main way that optimization algorithms (TopOpt) interact with a TopOptRep.  */
	virtual void getRealRep(std::vector<double>& realVec) const = 0;
	//! Copy the discrete topology representation into `realVec`
	virtual void getDiscreteRep(std::vector<int>& discVec) const = 0;
	//! Copies all optimization parameters into the input arguments.
	/*! The values returned in the arguments must fully specify the topology.  This function is used to
	 *  send values over MPI communication for parallel computations.
	 */
	virtual void getMPIRep(std::vector<std::vector<int> >& discreteVars, std::vector<std::vector<double> >& realVars) const = 0;
	//! Return the total number of optimization parameters
	virtual std::size_t getDataSize() const = 0;
	//! Return the number of optimization parameters in `sizes` if those parameters can be specified along more than 1 dimension
	virtual void getDataSize(std::vector<std::size_t>& sizes) const = 0;
	//! Return the spatial dimension of the TopOptRep (2 or 3)
	virtual unsigned getDimension() const = 0;
	//! Compute and return the volume fraction of the TopOptRep, useful for computing constraints.
	virtual double computeVolumeFraction() const = 0;
	//! Compute and return the gradient of the volume fraction of the TopOptRep
	virtual std::vector<double> computeGradVolumeFraction() const = 0;
	//! Computes the partial derivatives of the mesh element values with respect to the design variables
	/*! Returns a vector of maps in a sparse matrix-like representation.  Rows correspond to design variables and
   *  columns correspond to mesh element values.  
   */
	virtual std::vector<std::map<std::size_t, double>> diffRep() const = 0;
	//! Return the absolute expected magnitude of the optimization parameters
	/*! This is currently not implemented, but it may be necessary for problems with vastly different scales of 
	 *  parameters.  The scale can effect how finite differences and/or meshes are computed.
	 */
	virtual double getDataMagnitude() const {return 1.;}
	//! Return the parameters that define the settings of a TopOptRep, usually a collection of inputs from the input file
	/*! This is needed for sending new information about a TopOptRep over MPI and will be used in implementation constructors. */
	virtual void getDefiningParameters(std::vector<std::vector<int> >& discreteParams, 
										std::vector<std::vector<double> >& realParams) const = 0;
	//! Fix all values in `realVec` to the range required by this TopOptRep (usually 0-1)
	virtual void boundsCheck(std::vector<double>& realVec) const = 0;
	//! Return the name of this TopOptRep (ie. pixel, voxel, etc.)
	virtual std::string getName() const = 0;
	//! Return the TORType of this TopOptRep
	TORType getType() const {return myTORT;}
	//@}
	//! @name Data modification functions
	//@{
	//! Filters the data contained in `valVec` with a filter of a given radius
	/*! This is typically used to filter the gradient or densities in a pixel/voxel representation */
	virtual void filterData(std::vector<double>& valVec, double radius) const = 0;
	//! Filters the data contained in the TopOptRep implementation with a filter of a given radius
	virtual void filterData(double radius) = 0;
	//@}
protected:
	TORType myTORT;
};
}// namespace

#endif
