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

#include "genericmaterial.h"
#include "helper.h"

namespace Topologies{
using std::vector;
using std::ostream;

GenericMaterial::GenericMaterial()
{
}

GenericMaterial::GenericMaterial(vector<Real> inMat):
	constitutiveParams(inMat)
{
}

GenericMaterial::GenericMaterial(vector<Real> inMat, vector<Real> inMin, vector<Real> inMax):
	constitutiveParams(inMat),
	paramMin(inMin),
	paramMax(inMax)
{
	for(std::size_t k = 0; k < constitutiveParams.size(); k++)
		constitutiveParams[k] = MAX(paramMin[k], MIN(paramMax[k], constitutiveParams[k]));
}

GenericMaterial::GenericMaterial(vector<Real> inMat, vector<Real> inMin, vector<Real> inMax, vector<Real> pc):
	constitutiveParams(inMat),
	paramMin(inMin),
	paramMax(inMax),
	rgbPrintColor(pc)
{
	for(std::size_t k = 0; k < constitutiveParams.size(); k++)
		constitutiveParams[k] = MAX(paramMin[k], MIN(paramMax[k], constitutiveParams[k]));
}

GenericMaterial::GenericMaterial(std::vector<Real> inMat, const GenericMaterial& inMatLims):
	constitutiveParams(inMat),
	paramMin(inMatLims.paramMin),
	paramMax(inMatLims.paramMax)
{
}

GenericMaterial::GenericMaterial(std::vector<Real> inMat, const GenericMaterial& inMatLims, std::vector<Real> rgb):
	constitutiveParams(inMat),
	paramMin(inMatLims.paramMin),
	paramMax(inMatLims.paramMax),
	rgbPrintColor(rgb)
{
}

void GenericMaterial::swap(GenericMaterial& arg2)
{
	constitutiveParams.swap(arg2.constitutiveParams);
	paramMin.swap(arg2.paramMin);
	paramMax.swap(arg2.paramMax);
	rgbPrintColor.swap(arg2.rgbPrintColor);
}

bool GenericMaterial::operator==(const GenericMaterial& inel) const
{
	bool eq = true;
	for(std::size_t k = 0; k < constitutiveParams.size(); k++)
		eq &= constitutiveParams[k] == inel.constitutiveParams[k];
	return eq;
}

bool GenericMaterial::operator!=(const GenericMaterial& inel) const
{
	return !(*this == inel);
}

void GenericMaterial::genRandomMat()
{
	constitutiveParams.clear();
	for(std::size_t k = 0; k < paramMin.size(); k++)
		constitutiveParams.push_back(HelperNS::RandomGen::instance().randRealInRange(paramMin[k], paramMax[k]));
}

std::vector<Real> GenericMaterial::getPrintColor() const
{
	if(rgbPrintColor.size() == 0)
	{
		Real r = 0., g = 0., b = 0.;
		std::size_t nm = constitutiveParams.size();
		if(nm >= 1 && hasLimits())
		{
			Real mMin = getParameterMin(0), mMax = getParameterMax(0),
					 param = getParameter(0);
			r = (param - mMin)/(mMax - mMin);
		}
		if(nm >= 2 && hasLimits())
		{
			Real mMin = getParameterMin(1), mMax = getParameterMax(1),
					 param = getParameter(1);
			g = (param - mMin)/(mMax - mMin);
		}
		if(nm >= 3 && hasLimits())
		{
			Real mMin = getParameterMin(2), mMax = getParameterMax(2),
					 param = getParameter(2);
			b = (param - mMin)/(mMax - mMin);
		}
		std::vector<Real> tmpRGBPC;
		tmpRGBPC.push_back(r);
		tmpRGBPC.push_back(g);
		tmpRGBPC.push_back(b);
		return tmpRGBPC;
	}
	return rgbPrintColor;
}

Real GenericMaterial::dist(GenericMaterial mat2) const
{
	assert(constitutiveParams.size() == mat2.size());
	Real sum = 0.;
	for(std::size_t k = 0; k < constitutiveParams.size(); k++)
	{
		Real diff = constitutiveParams[k] - mat2.constitutiveParams[k];
		sum += diff*diff;
	}
	return sqrt(sum);
}

// Stream operator
ostream& operator<<(ostream& theStream, const GenericMaterial& inGM)
{
	for(std::size_t k = 0; k < inGM.constitutiveParams.size(); k++)
		theStream << ' ' << inGM.constitutiveParams[k];
	return theStream;
}


void GenericMaterial::setScaledMaterialParam(std::size_t k, Real val)
{
	assert(k < paramMax.size());
	assert(k < paramMin.size());
	assert(k < constitutiveParams.size());
	if(val > 1.)
		val = 1.;
	else if(val < 0.)
		val = 0.;
	constitutiveParams[k] = (paramMax[k] - paramMin[k])*val + paramMin[k];
}
}//namespace
