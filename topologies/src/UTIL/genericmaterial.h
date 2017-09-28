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

#ifndef GENERICMATERIAL_H
#define GENERICMATERIAL_H

#include <vector>
#include <iostream>
#include <assert.h>
#include "topologiesdefs.h"

namespace Topologies{
//! A basic class for implementing generic material properties
/*! This class can be used for storing and manipulating generic material properties.  
 *  It is mostly a wrapper for a std::vector, though minimum and maximum limits can
 *  be placed on the properties, as well as an RGB triplet for plotting purposes.
 */
class GenericMaterial
{
public:
	GenericMaterial();
	//! Constructor taking a vector representing a set of material properties
	GenericMaterial(std::vector<Real> inMat);
	//! Constructor taking a vector representing a set of material properties, as well as minimum and maximum values for those properties
	GenericMaterial(std::vector<Real> inMat, std::vector<Real> inMin, std::vector<Real> inMax);
	//! Constructor taking a vector representing a set of material properties, as well as minimum and maximum values for those properties and an RGB triplet for plotting
	GenericMaterial(std::vector<Real> inMat, std::vector<Real> inMin, std::vector<Real> inMax, std::vector<Real> pc);
	//! Constructor taking a vector representing a set of material properties, as well as minimum and maximum values for those properties, copied from inMatLims
	GenericMaterial(std::vector<Real> inMat, const GenericMaterial& inMatLims);
	//! Constructor taking a vector representing a set of material properties, as well as minimum and maximum values for those properties (copied from inMatLims) and an RGB triplet
	GenericMaterial(std::vector<Real> inMat, const GenericMaterial& inMatLims, std::vector<Real> rgb);
	GenericMaterial(const GenericMaterial& inEl);
	void swap(GenericMaterial& arg2);
	~GenericMaterial();

	//! Equality test, this tests only the material properties, not color or limits
	bool operator==(const GenericMaterial& inMat) const;
	//! Inquality test, this tests only the material properties, not color or limits
	bool operator!=(const GenericMaterial& inMat) const;
	//! Output stream insertion to print material properties
	friend std::ostream& operator<<(std::ostream& theStream, const GenericMaterial& inGM);
	//! Returns the kth parameter
	Real getParameter(std::size_t k) const;
	//! Returns the kth parameter minimum
	Real getParameterMin(std::size_t k) const;
	//! Returns the kth parameter maximum
	Real getParameterMax(std::size_t k) const;
	//! Returns the number of parameters
	std::size_t getNumParameters() const;
	//! Returns whether or not max and min limits are set for this GenericMaterial
	bool hasLimits() const;
	//! Returns the RGB triplet representing a print color
	/*! If no RGB color is defined, but the min and max values are, a color based on the values 
	 *  of the material properties is generated.
	 */
	std::vector<Real> getPrintColor() const;
	//! Replaces this GenericMaterial object with a random material with properties generated from a uniform distribution between the min & max values
	void genRandomMat();
	//! Returns the Euclidean distance between two GenericMaterial objects
	Real dist(GenericMaterial mat2) const;
	//! Sets the kth material property to scaled value `val` between the max and min 
	void setScaledMaterialParam(std::size_t k, Real val);
private:
	std::vector<Real> constitutiveParams, paramMin, paramMax;
	std::vector<Real> rgbPrintColor;
};

inline
Real GenericMaterial::getParameter(std::size_t k) const
{
	assert(k < constitutiveParams.size());
	return constitutiveParams[k];
}

inline
Real GenericMaterial::getParameterMin(std::size_t k) const
{
	assert(k < paramMin.size());
	return paramMin[k];
}

inline
Real GenericMaterial::getParameterMax(std::size_t k) const
{
	assert(k < paramMax.size());
	return paramMax[k];
}

inline
std::size_t GenericMaterial::getNumParameters() const 
{
	return constitutiveParams.size();
}

inline
bool GenericMaterial::hasLimits() const
{
	return(paramMin.size() > 0);
}
}// namespace
#endif
