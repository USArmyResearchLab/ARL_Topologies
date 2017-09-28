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

#ifndef ELEMENT_H
#define ELEMENT_H

#include <memory>
#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "point2d.h"
#include "UTIL/genericmaterial.h"
#include "globaldefs.h"

template <typename T> class ElemEdge;
class Cell;

enum PatchType
{
	ptLinTri = 1, ptLinQuad
};

//! Element is an abstract base class defining two dimensional geometric objects
/*! The class takes a template parameter defining the Point object to use to define
 *  its nodes.  In other words, the template paramter indicates whether or not the
 *  Element object lives in 2d or 3d space.  
 *  This class is the 2D equivalent of Cell and contains all geometric, topological,
 *  and numerical functionality to implement a 2D finite element method.
 */
template <class T>
class Element
{
public:
	//! Constructor taking a vector of defining nodes and a material property
	Element(const std::vector<T*>& ptVec, const Topologies::GenericMaterial& inMat);
	virtual ~Element();

	//! @name Element matrix computation functions
	//@{
	//! Returns the position in real space of the given area coordinate (\xi, \eta)
	virtual T eval(double xi, double eta) const = 0;
	//! Computes and returns the linear elastic element matrix
	virtual Eigen::MatrixXd getElemMat(const ElasticProblemType inEPT) const = 0;
	//! Computes and returns the Laplacian element matrix for given material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(double matVal) const = 0;
	//! Computes and returns the Laplacian element matrix for given anisotropic material value matVal
	virtual Eigen::MatrixXd getLaplacianElemMat(const Eigen::MatrixXd& matVal) const = 0;
	//! Computes and returns the mass matrix (integral of basis * test functions)
	virtual Eigen::MatrixXd getMassMat() const = 0;
	//! Computes and returns the integral of the testing functions
	virtual Eigen::MatrixXd getTFInt() const = 0;
	//! Computes and returns the thermal expansion stiffness matrix
	virtual Eigen::MatrixXd getThermalExpansionMat(const ElasticProblemType inEPT, double alphaTE) const = 0;
	//! Returns a vector containing the basis function values at the given local area coordinate (\xi, \eta)
	virtual Eigen::VectorXd getInterpolationVec(double xi, double eta) const = 0;
	//! Returns a matrix containing the gradient of the basis functions at element k, computed at the centroid
	virtual Eigen::MatrixXd getGradBF() const = 0;
	//! Computes and returns the Laplacian of the basis function: div(grad(phi))
	/*! Note that these functions may be inaccurate and should be thoroughly tested before use.*/
	virtual Eigen::MatrixXd getLaplacianBF1() const = 0;
	//! Computes and returns the Laplacian of the basis function: grad(grad(phi))
	/*! Note that these functions may be inaccurate and should be thoroughly tested before use.*/
	virtual Eigen::MatrixXd getLaplacianBF2() const = 0;
	//@}
	//! @name Geometric functions
	//@{
	//! Returns the area of the Element
	virtual double getArea() const = 0;
	//! Returns the unitary basis vector associated with the first area coordinate (\xi) at the point (\xi, \eta).
	/*! Unitary basis vectors are the derivative of the position vector with respect to the area coordinates */
	virtual T getXiVector(double xi, double eta) const = 0;
	//! Returns the unitary basis vector associated with the second area coordinate (\eta) at the point (\xi, \eta).
	/*! Unitary basis vectors are the derivative of the position vector with respect to the area coordinates */
	virtual T getEtaVector(double xi, double eta) const = 0;
	//! Returns the Jacobian of the element at the point (\xi, \eta)
	virtual double getJacobian(double xi, double eta) const = 0;
	//! Returns the centroid of the Element, or the average of the nodal coordinates
	T centroid() const;
	//@}
	//! @name Topological functions
	//@{
	//! Returns the sign of the Jacobian
	virtual int getOrientation() const = 0;
	//! Swaps the defining area coordinates if necessary.  Ensure positivity of the Jacobian
	virtual void switchXiEta();
	//! Test for equality of two Elements.  Elements are equal if they contain the same set of nodes, even in any order.
	bool operator==(const Element& arg2) const;
	//! Returns an the kth edge of the Element
	ElemEdge<T>* getEdge(unsigned k) const;
	//! Returns an the kth node of the Element
	T* getNode(unsigned k) const;
	//! Returns the number of nodes associated with the Element (eg. 3 for a triangle)
	unsigned getNumNodes() const {return nodeVec.size();}
	//! Sets the argument ptNew as the kth node
	void setNode(T* ptNew, unsigned k);
	//! Creates a a new edge with nodes determined by the argument edgeNum
	std::unique_ptr<ElemEdge<T>> createEdge(unsigned short int edgeNum) const;
	//! Adds a new edge to the Element
	void addEdge(ElemEdge<T>* theEdge);
	//! Returns whether or not this Element's nodes have been reordered to set a standard edge/node ordering 
	bool isReordered() const;
	//! Sets a flag indicating that this Element's nodes should not be reordered again
	void hasBeenReordered();
	//! Reorders the edges and nodes of an element so that they are in a standard order.
	/*! For example, a triangle has a natural ordering such that the 0th edge is opposite the 0th node.  
	 *  Non-simplexes can have other specified orders
	 */
	void setStandardEdgeOrder();
	//! Adds nodal basis function id bfGID to this element, associated with node theNode
	void addNodeBases(const T* theNode, std::size_t bfGID);
	//! Returns a boundary edge associated with this element, if one exists
	ElemEdge<T>* findBoundaryEdge() const;
	//@}
	//! @name Data access
	//@{
	//! Returns an Element's material properties
	Topologies::GenericMaterial getGenericMaterial() const;
	//! Returns the Element's type
	PatchType getPatchType() const {return myPT;}
	//! Returns whether or not a point is contained in this Element, for both 2d and 3d.
	/*! Note that this is a topological test, not a geometric test, i.e the node's pointer is checked, not its coordinates.
	 */
	bool containsNode(const T* tstNode) const;
	//! Returns the global basis function number for local basis function locBF
	std::size_t getGlobalBF(unsigned locBF) const;
	//! Returns the global basis function for a local basis function number that exists on a boundary
	std::size_t getGlobalBFonBoundary(unsigned locBF) const;
	//! Returns the edge number for a given ElemEdge
	unsigned getLocEdgeID(ElemEdge<T>* findEdge) const;
	//@}
	//! @name Functions used with the Cell class
	//@{
	//! Associates Cell inCell with this Element, and is used to define a Cell's face
	bool addCell(Cell* inCell);
	//! Returns whether or not this Element is on the boundary of a volume
	bool isBoundaryPatch() const {return pCell1 == nullptr;}
	//! Returns a Cell.  This will always be a well-defined Cell
	Cell* cell0() const {return pCell0;}
	//! Returns a Cell.  This may be a nullptr if this Element is on the boundary
	Cell* cell1() const {return pCell1;}
	//@}
protected:
	std::pair<double,double> getCijs(const ElasticProblemType inEPT) const;

	PatchType myPT;
	Topologies::GenericMaterial material;
	std::vector<T*> nodeVec;
	std::vector<ElemEdge<T>*> edgeVec;
	std::vector<std::size_t> bfVec;
	Cell *pCell0, *pCell1;
	bool reordered;

	static const unsigned cRho = 0, cLambda = 1, cMu = 2;
};

//! Hash function for testing Element equality, used in std::unordered_set
/*! Adapted from http://stackoverflow.com/questions/1536393/good-hash-function-for-permutations 
 *  For this class, we want to ensure that elements that contain permutations of nodes test
 *  as equal, so a hash that gives the same value for any permutation of node pointers is needed.
 */
template<typename T>
struct Element_hash
{
	typedef const Element<T>* const argument_type;
	typedef std::size_t result_type;
	result_type operator()(argument_type const& s) const
	{
		const std::size_t R = 1779033703;
		std::size_t h1 = 1.;
		for(unsigned k = 0; k < s->getNumNodes(); ++k)
			h1 *= (R + 2*reinterpret_cast<size_t>(s->getNode(k)));
		return h1/2;
	}
};

//! Equality test for elements, needed for using std::unordered_set
template<typename T>
struct Element_equal
{
	typedef const Element<T>* const argument_type;
	bool operator()(argument_type const& x, argument_type const& y) const
	{
		return *x == *y;
	}
};

template <typename T> 
ElemEdge<T>* Element<T>::getEdge(unsigned k) const
{
	assert(k < edgeVec.size());
	return edgeVec[k];
}

template <typename T> 
T* Element<T>::getNode(unsigned k) const
{
	return nodeVec[k];
}

template <typename T>
void Element<T>::setNode(T* ptNew, unsigned k)
{
	nodeVec[k] = ptNew;
}

template <typename T> 
bool Element<T>::isReordered() const 
{
	return reordered; 
}

template <typename T> 
void Element<T>::hasBeenReordered()
{ 
	reordered = true; 
}

template <typename T>
std::size_t Element<T>::getGlobalBF(unsigned locBF) const
{
	assert(locBF < bfVec.size());
	return bfVec[locBF];
}

template <typename T>
Topologies::GenericMaterial Element<T>::getGenericMaterial() const
{
	return material;
}

#endif
