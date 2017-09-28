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

#ifndef TOMESHPROCESSING_H
#define TOMESHPROCESSING_H

#include "cgal_types.h"

namespace Topologies{
class TOMesh;

//! A collection of functions for processing of TOMesh objects
/*! Processing includes computing element areas/volumes and centroids. 
*/
namespace TOMeshProcessing
{
	//! Returns the volume of the element with id `elemID` in `inMesh`
	double computeElementVolume(std::size_t elemID, const TOMesh* const inMesh);
	//! Returns the centroid of the element with id `elemID` in `inMesh`, assuming a 2d mesh
	Point_2_base getElementCentroid2D(std::size_t elemID, const TOMesh* const inMesh);
	//! Returns the centroid of the element with id `elemID` in `inMesh`, assuming a 3d mesh
	Point_3_base getElementCentroid3D(std::size_t elemID, const TOMesh* const inMesh);
}
}
#endif

