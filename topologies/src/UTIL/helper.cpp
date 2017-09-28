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

#include "helper.h"

namespace Topologies{
int HelperNS::round(Real x)
{
	return x > 0. ? static_cast< int >(x + 0.5) : static_cast< int >(x - 0.5);
}

Uint HelperNS::roundToUint(Real x)
{
	assert(x >= 0.);
	return static_cast<Uint>(x + 0.5);
}

std::string HelperNS::upperCase(std::string instring)
{
	std::transform(instring.begin(), instring.end(), instring.begin(), toupper);
	return instring;
}

std::string HelperNS::lowerCase(std::string instring)
{
	std::transform(instring.begin(), instring.end(), instring.begin(), tolower);
	return instring;
}

std::vector<double> HelperNS::getPenalizedPixels(const std::vector<double>& pixelArray, double penalPower, double minDensity)
{
	HelperNS::powPenal powFunc(penalPower);
	std::vector<double> penalPixels(pixelArray.size());
	std::transform(pixelArray.begin(), pixelArray.end(), penalPixels.begin(), powFunc);
	HelperNS::affineTrans1d at1d((1. - minDensity), minDensity);
	std::transform(penalPixels.begin(), penalPixels.end(), penalPixels.begin(), at1d);
	return penalPixels;
}
}

