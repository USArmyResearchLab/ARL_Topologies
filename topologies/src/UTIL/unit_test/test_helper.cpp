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

#define CATCH_CONFIG_MAIN

#include "helper.h"
#include "catch.hpp"
#include <iostream>
#include <string>

using namespace Topologies;

TEST_CASE("Testing functions in HelperNS namespace","[HelperNS]")
{
	using namespace HelperNS;
	// Random functions
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randRealInRange<double>() <= 1.);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randRealInRange<double>() >= 0.);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randRealInRange<double>(-1., 0.) <= 0.);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randRealInRange<double>(-1., 0.) >= -1.);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randIntInRange<int>(-3, -1) <= -1);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randIntInRange<int>(-3, -1) >= -3);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randIntInRange<unsigned>(1, 2) <= 2);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randIntInRange<unsigned>(1, 2) >= 1);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randIntInRange<std::size_t>(100, 2000) <= 2000);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(RandomGen::instance().randIntInRange<std::size_t>(100, 2000) >= 100);
	// Round
	double dx = 0.1;
	for(double val = -1.; val <= -0.5; val += 0.1)
		REQUIRE(HelperNS::round(val) == -1);
	for(double val = 0.5; val <= 0.; val += 0.1)
		REQUIRE(HelperNS::round(val) == 0);
	for(double val = 0.; val < 0.5; val += 0.1)
		REQUIRE(HelperNS::round(val) == 0);
	for(double val = 0.5; val <= 1.; val += 0.1)
		REQUIRE(HelperNS::round(val) == 1);
	for(double val = 0.; val < 0.5; val += 0.1)
    REQUIRE(HelperNS::roundToUint(val) == 0);
  for(double val = 0.5; val <= 1.; val += 0.1)
    REQUIRE(HelperNS::roundToUint(val) == 1);
	// Text processing
	std::string testStr("This Is A TesT STring...");
	REQUIRE(upperCase(testStr) == "THIS IS A TEST STRING...");
	REQUIRE(lowerCase(testStr) == "this is a test string...");
	// Penalized pixels
	std::vector<double> array1 = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
	double penalPower = 2.5;
	double minDensity = 0.1;
	std::vector<double> res = getPenalizedPixels(array1, penalPower, minDensity);
	REQUIRE(res.size() == array1.size());
	for(std::size_t k = 0; k < array1.size(); ++k)
		REQUIRE(res[k] == Approx((1. - minDensity)*pow(array1[k], penalPower) + minDensity));
	// RGWrapper
	RGWrapper testRGW(std::make_pair(-1., 1.));
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(testRGW() <= 1.);
	for(unsigned k = 0; k < 100; ++k)
		REQUIRE(testRGW() >= -1.);
	// r2d
	r2d testr2d(0.5);
	for(double val = 0.; val < 0.5; val += 0.1)
		REQUIRE(testr2d(val) == 0);
	for(double val = 0.5; val <= 1.; val += 0.1)
		REQUIRE(testr2d(val) == 1);
	// affineTrans1d
	affineTrans1d testat1d(2., 3.);
	REQUIRE(testat1d(1.) == Approx(5.));
	REQUIRE(testat1d(-1.) == Approx(1.));
	// powPenal
	powPenal testpowPenal(2.);
	REQUIRE(testpowPenal(3.3) == Approx(3.3*3.3));
	powPenal dtestpowPenal(2., true);
	REQUIRE(dtestpowPenal(3.3) == Approx(2.*3.3));
	// powPenalMin
	powPenalMin testpowPenalMin(3., 0.1);
	REQUIRE(testpowPenalMin(3.3) == Approx(3.3*3.3*3.3*0.9 + 0.1));
	REQUIRE(testpowPenalMin(0.) == Approx(0.1));
	powPenalMin dtestpowPenalMin(3., 0.1, true);
	REQUIRE(dtestpowPenalMin(3.3) == Approx(3.*3.3*3.3*0.9));
	REQUIRE(dtestpowPenalMin(0.) == Approx(0.));
	// powSum
	powSum testpowSum(3.);
	REQUIRE(testpowSum(4., 5.) == Approx(4. + 5.*5.*5.));
	// regularizedHeaviside
	regularizedHeaviside testrh(9.);
	REQUIRE(testrh(0.5) == Approx(9.889527083638011e-01));
	regularizedHeaviside dtestrh(9., true);
	REQUIRE(dtestrh(0.5) == Approx(1.001043786482674e-01));
	// thresholdHeaviside
	thresholdHeaviside testth(0.25, 11.);
	REQUIRE(testth(0.5) == Approx(9.959132968164297e-01));
	thresholdHeaviside dtestth(0.25, 11., true);
	REQUIRE(dtestth(0.5) == Approx(8.954303588654487e-02));
	// divVal
	divVal testdv(2.);
	REQUIRE(testdv(1.) == Approx(0.5));
	// multVal
	multVal<int> testmv(2);
	REQUIRE(testmv(1) == 2);
	// greaterThanX
	greaterThanX testgtx(0.5);
	REQUIRE(testgtx(1.));
	REQUIRE(!testgtx(0.));
	// lessThanX
	lessThanX testltx(0.5);
	REQUIRE(!testltx(1.));
	REQUIRE(testltx(0.));
	// vector norm
	REQUIRE(norm(array1) == Approx(1.962141687034859));
	// vecSum
	std::vector<double> array2 = {0., -0.1, 0.2, -0.3, 0.4, -0.5, 0.6, -0.7, 0.8, -0.9, 1.};
	std::vector<double> exactRes = {0., 0., 0.4, 0., 0.8, 0., 1.2, 0., 1.6, 0., 2.};
	res = vecSum(array1, array2);
	REQUIRE(res.size() == array1.size());
	REQUIRE(res.size() == array2.size());
	REQUIRE(res.size() == exactRes.size());
	for(std::size_t k = 0; k < res.size(); ++k)
		REQUIRE(res[k] == Approx(exactRes[k]));
	// VecMinus
	exactRes = {0., 0.2, 0., 0.6, 0., 1., 0., 1.4, 0., 1.8, 0.};
	res = vecMinus(array1, array2);
	REQUIRE(res.size() == array1.size());
	REQUIRE(res.size() == array2.size());
	REQUIRE(res.size() == exactRes.size());
	for(std::size_t k = 0; k < res.size(); ++k)
		REQUIRE(res[k] == Approx(exactRes[k]));
	// vecScalarMult
	exactRes = {0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.};
	res = vecScalarMult(2., array1);
	REQUIRE(res.size() == array1.size());
	REQUIRE(res.size() == array2.size());
	REQUIRE(res.size() == exactRes.size());
	for(std::size_t k = 0; k < res.size(); ++k)
		REQUIRE(res[k] == Approx(exactRes[k]));
}
