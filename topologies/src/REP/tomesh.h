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

#ifndef TOMESH_H
#define TOMESH_H

#include <vector>
#include <map>
#include "cgal_types.h"

namespace Topologies{
//! Simple polymorphic wrapper struct to hold a mesh and give access to its elements, nodes, and material properties
class TOMesh
{
public:
	//! Constructor taking the dimension of the mesh (2 or 3)
	TOMesh(unsigned inDimNum) : m_dimNum(inDimNum) {}
	TOMesh(unsigned inDimNum, std::vector<std::vector<std::size_t>> inCVec);
	TOMesh(unsigned inDimNum, std::vector<std::vector<std::size_t>> inCVec, std::vector<double> inOVec);
	TOMesh(unsigned inDimNum, std::vector<std::vector<std::size_t>> inCVec, std::vector<int> inMIDVec);
	TOMesh(unsigned inDimNum, std::vector<std::vector<std::size_t>> inCVec, std::vector<double> inOVec, 
		std::vector<int> inMIDVec);
	virtual ~TOMesh() {}
	virtual std::unique_ptr<TOMesh> clone() const = 0;

	//! Mesh dimension
	unsigned dimNum() const {return m_dimNum;}
	//! Returns the number of elements in the mesh
	std::size_t getNumElements() const {return connectivityVec.size();}
	//! Returns the number of nodes (points) in the mesh
	virtual std::size_t getNumNodes() const = 0;
	//! Returns the element connectivity (set of defining node IDs) for element `kelem` 
	const std::vector<std::size_t>& getElementConnectivity(std::size_t kelem) const {return connectivityVec[kelem];}
	//! Returns the scalar optimization value for element kelem
	/*! This value is typically used to scale material properties. */
	double getOptVal(std::size_t kelem) const {return optVals[kelem];}
	//! Returns the material ID for `kelem`
	int getMatID(std::size_t kelem) const {return matIDVec[kelem];}
	//! Returns the node (point) with ID `kn` (for 2d meshes)
	virtual Point_2_base getNode2D(std::size_t kn) const = 0;
	//! Returns the node (point) with ID `kn` (for 3d meshes)
	virtual Point_3_base getNode3D(std::size_t kn) const = 0;
	//! Sets the optimization values for the mesh (densities)
	void setOptVals(const std::vector<double>& inVec);
	//! Sets one optimization value
	void setOptVal(std::size_t ke, double val);
	//! Sets the material ids for the mesh, which may be useful for a FEM solver
	void setMatIDs(const std::vector<int>& inVec);
	//! Sets one material id for the mesh
	void setMatIDs(std::size_t ke, int val);
protected:
	//! Swaps this TOMesh with arg
	void swap(TOMesh& arg);
	//! A vector containing all element connectivity in the mesh
	std::vector<std::vector<std::size_t>> connectivityVec;
	//! A vector containing all optimization values
	std::vector<double> optVals;
	//! A vector containing all material ids, used for specifiying materials in objective functions
	std::vector<int> matIDVec;
	//! Mesh dimension
	unsigned m_dimNum;
};

//! A 2d mesh, derived from TOMesh
class TOMesh2D : public TOMesh
{
public:
	TOMesh2D() : TOMesh(2) {}
	//! Constructor that converts a CGAL Constrained Delaunay Triangulation to a TOMesh2D
	TOMesh2D(const CDT_2& inMesh);
	//! Constructor that copies or moves the arguments `inPtVec` and `inConnVec` to construct a TOMesh2D
	TOMesh2D(std::vector<Point_2_base> inPtVec, std::vector<std::vector<std::size_t>> inConnVec);
	//! Constructor that copies or moves the arguments `inPtVec`, `inConnVec`, and `inOptVec` to construct a TOMesh
	TOMesh2D(std::vector<Point_2_base> inPtVec, std::vector<std::vector<std::size_t>> inConnVec, std::vector<double> inOptVec);
	//! Constructor that copies or moves the arguments `inPtVec`, `inConnVec`, and `inMIDVec` to construct a TOMesh
	TOMesh2D(std::vector<Point_2_base> inPtVec, std::vector<std::vector<std::size_t>> inConnVec, std::vector<int> inMIDVec);
	//! Constructor that copies or moves the arguments `inPtVec`, `inConnVec`, `inOptVec`, and `inMIDVec` to construct a TOMesh
	TOMesh2D(std::vector<Point_2_base> inPtVec, std::vector<std::vector<std::size_t>> inConnVec,
		std::vector<double> inOptVec, std::vector<int> inMIDVec);
	virtual ~TOMesh2D() {}
	//! Swaps this TOMesh2D wtih arg
	void swap(TOMesh2D& arg);
	virtual std::unique_ptr<TOMesh> clone() const	{return std::unique_ptr<TOMesh>(new TOMesh2D(*this));	}
	
	virtual std::size_t getNumNodes() const {return ptVec.size();}
	virtual Point_2_base getNode2D(std::size_t kn) const	{return ptVec[kn];}
	virtual Point_3_base getNode3D(std::size_t kn) const {return Point_3_base(0., 0., 0.);} // Not implemented
private:
	void setupFromCDT2(const CDT_2& inMesh);
	//! A vector containing all points (nodes) in the mesh
  std::vector<Point_2_base> ptVec;
};

//! A 3d mesh, derived from TOMesh
class TOMesh3D : public TOMesh
{
public:
	TOMesh3D() : TOMesh(3) {}
	//! Constructor that copies or moves the arguments `inPts` and `inElems` to construct a TOMesh3D
	/* Note that edges and faces are not set up with this constructor! */
	TOMesh3D(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems);
	//! Constructor that copies or moves the arguments `inPts`, `inElems`, and `inOptVals` to construct a TOMesh3D
	/*! Note that edges and faces are not set up with this constructor! */
	TOMesh3D(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems, std::vector<double> inOptVals);
	//! Constructor that copies or moves the arguments `inPts`, `inElems`, and `inMIDVec` to construct a TOMesh3D
	/*! Note that edges and faces are not set up with this constructor! */
	TOMesh3D(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems, std::vector<int> inMIDVec);
	//! Constructor that copies the arguments `inPts`, `inElems`, `inOptVec`, and `inMIDVec` to construct a TOMesh3D
	/*! Note that edges and faces are not set up with this constructor! */
	TOMesh3D(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems,
		std::vector<double> inOptVals, std::vector<int> inMIDVec);
	//! Constructor that converts a CGAL 3d tetrahedral mesh to TOMesh3D
	TOMesh3D(const C3t3& inMesh);
	virtual ~TOMesh3D() {}
	//! Swaps this TOMesh3D wtih arg
	void swap(TOMesh3D& arg);
	virtual std::unique_ptr<TOMesh> clone() const {return std::unique_ptr<TOMesh>(new TOMesh3D(*this)); }
	
	//! A vector containing the mesh edges in terms of node connectivity
	std::vector<std::vector<std::size_t> > edgeVec;
	//! A vector containing the mesh faces in terms of node connectivity
	std::vector<std::vector<std::size_t> > faceVec;

	virtual std::size_t getNumNodes() const {return ptVec.size();}
	virtual Point_2_base getNode2D(std::size_t kn) const {return Point_2_base(0., 0.);} // Not implemented
	virtual Point_3_base getNode3D(std::size_t kn) const { return ptVec[kn];}
private:
	void setupFromC3t3(const C3t3& inMesh);

	std::vector<Point_3_base> ptVec;
};

//! A 3d surface mesh, derived from TOMesh
class TOMesh3DSurface : public TOMesh
{
private:
	typedef typename C2t3::Triangulation TOM_Tr;
	typedef typename TOM_Tr::Vertex_handle TOM_Vertex_handle;
public:
	//! Constructor that converts a CGAL C2t3 object to a TOMesh
	TOMesh3DSurface(const C2t3& inMesh);
	//! Constructor that copies or moves a set of points and connectivity
	TOMesh3DSurface(std::vector<Point_3_base> inPts, std::vector<std::vector<std::size_t>> inElems);
	virtual ~TOMesh3DSurface() {}
  //! Swaps this TOMesh3DSurface wtih arg
  void swap(TOMesh3DSurface& arg);
	virtual std::unique_ptr<TOMesh> clone() const {return std::unique_ptr<TOMesh>(new TOMesh3DSurface(*this)); }	

	virtual std::size_t getNumNodes() const {return ptVec.size();}
	virtual Point_2_base getNode2D(std::size_t kn) const {return Point_2_base(0., 0.);} // Not implemented
	virtual Point_3_base getNode3D(std::size_t kn) const {return ptVec[kn];}
private:
	std::vector<Point_3_base> ptVec;

	void setupFromC2t3(const C2t3& inMesh);
};
} // namespace
#endif

