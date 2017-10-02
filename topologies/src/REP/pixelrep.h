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

#ifndef PIXELREP_H
#define PIXELREP_H

#include "volmesh2d.h"
#include "inputloaderrep.h"
#include <memory>

namespace Topologies{
class GenericMaterial;

//! A topology representation that discretizes a rectangular region into pixels
/*! This is the most basic topology representation.  A region is broken into rectangular pixels
 *  each of which has an optimization parameter.  Note that in order for this to be effective
 *  sensitivity filtering is necessary.  Meshes can use either tri or quad elements.
 *  This class further specializes VolMesh2D (which can be used for arbitrary 2d meshes).  It is also
 *  parameterized by a penalization function and a projection function.
 */
template <typename PenaltyFunc = HelperNS::powPenalMin, typename ProjectionFunc = HelperNS::defaultProjFunc>
class PixelRep : public VolMesh2D<PenaltyFunc, ProjectionFunc>
{
public:
	PixelRep(TORType inTORT, const InputLoader::TORGenericVolume& inputParams, const std::vector<double>& inPenalParams,
    const std::vector<double>& inProjParams);
	PixelRep(TORType inTORT, const std::vector<std::vector<int> >& discreteParams, 
			const std::vector<std::vector<double> >& realParams);
	virtual ~PixelRep();
	PixelRep(const PixelRep& copy);
	PixelRep(PixelRep&& copy) : VolMesh2D<PenaltyFunc, ProjectionFunc>(copy.myTORT) {swap(copy);}
	PixelRep& operator=(PixelRep rhs){swap(rhs); return *this;}
	void swap(PixelRep& arg2);
	virtual std::unique_ptr<TopOptRep> clone() const;

	//! @name Structural modification functions
	//@{
	virtual void refine();
	virtual void prune();
	//@}
	//! @name Data access
	//@{
	virtual void getDataSize(std::vector<std::size_t>& sizes) const;
	virtual void getDefiningParameters(std::vector<std::vector<int> >& discreteParams,
                                     std::vector<std::vector<double> >& realParams) const;
	virtual std::string getName() const {return getClassName();}
	//! Returns TopOptRep implementation's class name: pixel
	static std::string getClassName() {return "pixel";}
	//@}

private:
	void finishSetup();
	void refineElementUnknowns();
	void refineNodalUnknowns();

	typedef VolMesh<PenaltyFunc, ProjectionFunc> VM;
	typedef VolMesh2D<PenaltyFunc, ProjectionFunc> VM2D;

	double width, height; // Region physical size
	unsigned nx, ny;
	MeshElementType myMET;
};

template <typename PenaltyFunc, typename ProjectionFunc>
PixelRep<PenaltyFunc, ProjectionFunc>::PixelRep(const PixelRep<PenaltyFunc, ProjectionFunc>& copy) :
	VolMesh2D<PenaltyFunc, ProjectionFunc>(copy),
	nx(copy.nx),
	ny(copy.ny),
	width(copy.width),
	height(copy.height),
	myMET(copy.myMET)
{
}

template <typename PenaltyFunc, typename ProjectionFunc>
void PixelRep<PenaltyFunc, ProjectionFunc>::swap(PixelRep<PenaltyFunc, ProjectionFunc>& arg2)
{
	VM2D::swap(arg2);
	std::swap(nx, arg2.nx);
	std::swap(ny, arg2.ny);
	std::swap(width, arg2.width);
	std::swap(height, arg2.height);
	std::swap(myMET, arg2.myMET);
}
}
#endif

