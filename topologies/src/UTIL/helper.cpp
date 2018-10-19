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
#include <fstream>

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

bool HelperNS::isFileReadable(std::string const& fileName)
{
	std::ifstream ifs(fileName);
	return ifs.good();
}

std::vector<double> HelperNS::SparseMatrix::operator*(std::vector<double> const& vec) const
{
	std::vector<double> res(m_mat.size(), 0.);
	auto rowIt = res.begin();
	for(auto const& rowVec : m_mat)
	{
		for(auto const& colPair : rowVec)
			*rowIt += value(colPair)*vec[index(colPair)];
		++rowIt;
	}
	return res;
}

double HelperNS::SparseMatrix::operator()(std::size_t row, std::size_t col) const
{
	SparseRow const& matRow = m_mat[row];
	for(auto const& colPair : matRow)
	{
		if(index(colPair) == col)
			return value(colPair);
	}
	return 0.;
}

std::vector<double> HelperNS::SparseMatrix::transposeTimes(std::vector<double> const& vec, std::size_t ncols) const
{
	std::vector<double> res(ncols, 0.);
	for(std::size_t krow = 0; krow < m_mat.size(); ++krow)
		for(auto const& colPair : m_mat[krow])
			res[index(colPair)] += vec[krow]*value(colPair);
	return res;
}

}

