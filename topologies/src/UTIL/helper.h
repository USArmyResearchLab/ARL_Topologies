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

#ifndef HELPER_NS_H
#define HELPER_NS_H

#include <algorithm>
#include <functional>
#include <numeric>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <random>
#include "topologiesdefs.h"

namespace Topologies{
//! A collection of utilities (functions and structs) for aiding more important things
/*! This namespace is a set of functions and structs to help other classes and functions.
*  Some functionality includes some basic structs to use with std::transform for applying
*  functions to a container of scalars and some wrappers for randomness.
*/
namespace HelperNS
{
	//! Singleton wrapper class for C++ random functionality
	class RandomGen
	{
	public:
		//! Returns the only instance of RandomGen
		static RandomGen& instance()
		{
			static RandomGen smInstance;
			return smInstance;
		}
		//! Resets the seed of the random number engine
		void reseed()
		{
			std::random_device r;
			dre = std::default_random_engine(r());
		}
		//! Sets the seed of the random number engine to specific value
		void reseed(double seed)
		{
			dre = std::default_random_engine(seed);
		}
		//! Generates a random integer from a uniform distribution between given limits
		template <typename T>
		T randIntInRange(T left, T right)
		{
			std::uniform_int_distribution<T> uniform_dist(left, right);
			return uniform_dist(dre);
		}
		//! Generates a random real from a uniform distribution between two limits
		//! Default range is 0 and 1
		template <typename T>
		T randRealInRange(T left = 0., T right = 1.)
		{
			std::uniform_real_distribution<T> uniform_dist(left, right);
			return uniform_dist(dre);
		}
		//! Returns the result of weighted coin flip
		bool coinFlip(double prob = 0.5)
		{
			return randRealInRange<double>() < prob;
		}
	private:
		RandomGen()
		{
			std::random_device r;
			dre = std::default_random_engine(r());
		}
		RandomGen(const RandomGen& copy){}
		RandomGen& operator=(RandomGen arg2){return *this;}
	private:
		std::default_random_engine dre;
	};
	//! Indirection wrapper for RandomGen, so that it can be used with std::generate, etc.
	class RGWrapper
	{
	public:
		RGWrapper() : randRange(std::make_pair(0.,1.)) {}
		RGWrapper(const std::pair<double,double>& inRange) : randRange(inRange) {}
		//! Returns a random real number within randRange
		double operator()() {return RandomGen::instance().randRealInRange(randRange.first, randRange.second);}
	private:
		std::pair<double, double> randRange;
	};
	//! Round to the nearest integer.
	int round(Real x);
	//! Round to the nearest unsigned integer.
	Uint roundToUint(Real x);
	//! Transform argument STL string to uppercase
	std::string upperCase(std::string instring);
	//! Transform argument STL string to lowercase
	std::string lowerCase(std::string instring);
	//! Applies a power law penalization to the input vector
	std::vector<double> getPenalizedPixels(const std::vector<double>& pixelArray, double penalPower, double minDensity);
	//! Struct to convert reals between 0 and 1 to ints with a threshold value
	struct r2d
	{
		double threshold;
		explicit r2d(double inThreshold) : threshold(inThreshold){}
		int operator()(double inval){return inval < threshold ? 0 : 1;}
	};
	//! Struct to apply an affine (a*x + b) transformation 
	struct affineTrans1d
	{
		double a, b;
		affineTrans1d(double ina, double inb) : a(ina), b(inb){}
		double operator()(double x){return a*x + b;}
	};
	//! Struct to apply a power law penalization
	struct powPenal
	{
		const double penal;
		const bool diff;
		explicit powPenal(double inPenal, bool inDiff = false) : penal(inPenal), diff(inDiff){}
		powPenal(const std::vector<double>& itl, bool inDiff = false) : penal(*itl.begin()), diff(inDiff) {assert(itl.size() == 1);}
		double operator()(double x){return diff ? g(x) : f(x);}
		double f(double inval){return pow(inval, penal);}
		double g(double inval){return penal*pow(inval, penal-1.);}
	};
	//! Struct to apply a power law with a minimum value
	struct powPenalMin
	{
		const double penal, minVal;
		const bool diff;
		powPenalMin(double inPenal, double inMinVal, bool inDiff = false) : penal(inPenal), minVal(inMinVal), diff(inDiff) {}
		powPenalMin(const std::vector<double>& itl, bool inDiff = false) : penal(*itl.begin()), minVal(*(itl.begin()+1)), diff(inDiff) 
		{assert(itl.size() == 2);}
		double operator()(double x){return diff ? g(x) : f(x);}
		double f(double x){return pow(x, penal)*(1 - minVal) + minVal;}
		double g(double x){return penal*pow(x, penal - 1.)*(1 - minVal);}
	};
	//! Struct to apply pow and sum
	struct powSum
	{
		double exponent;
		explicit powSum(double inExp) : exponent(inExp){}
		double operator()(const double& a, const double& b) {return a + pow(b, exponent);}
	};
	//! Struct to apply a regularized Heaviside function
	struct regularizedHeaviside
	{
		const double beta;
		const bool diff;
		explicit regularizedHeaviside(double inBeta, bool inDiff = false) : beta(inBeta), diff(inDiff) {}
		regularizedHeaviside(const std::vector<double>& itl, bool inDiff = false) : beta(*itl.begin()), diff(inDiff) 
		{assert(itl.size() == 1);}
		double operator()(double x){return diff ? g(x) : f(x);}
		double f(double inval){return 1. - exp(-beta*inval) + inval*exp(-beta);}
		double g(double inval){return beta*exp(-beta*inval) + exp(-beta);}
	};
	//! Struct to apply a shifted regularized Heaviside function
	struct thresholdHeaviside
	{
		const double T, beta;
		const bool diff;
		thresholdHeaviside(double inT, double inBeta, bool inDiff = false) : T(inT), beta(inBeta), diff(inDiff) {}
		thresholdHeaviside(const std::vector<double>& itl, bool inDiff = false) : T(*itl.begin()), beta(*(itl.begin()+1)), diff(inDiff) 
		{assert(itl.size() == 2);}
		double operator()(double x){return diff ? g(x) : f(x);}
		double f(double x){return (tanh(beta*T) + tanh(beta*(x - T)))/(tanh(beta*T) + tanh(beta*(1. - T)));}
		double g(double x){return (1 - tanh(beta*(x - T))*tanh(beta*(x - T)))*beta/(tanh(beta*T) + tanh(beta*(1. - T)));}
	};
	//! Struct to do nothing to a value, used as a default template parameter for VolMesh2D
	struct defaultProjFunc
	{
		const bool diff;
		explicit defaultProjFunc(bool inDiff = false) : diff(inDiff) {}
		defaultProjFunc(const std::vector<double>& itl, bool inDiff = false) : diff(inDiff) {}
		double operator()(double x){return diff ? g(x) : f(x);}
		double f(double x){return x;}
		double g(double x){return 1.;}
	};
	//! Struct for a linear hat function, used in Filter2D and Filter3D
	struct linearHat
	{
		const double rad;
		explicit linearHat(double inRad) : rad(inRad) {}
		double operator()(double x)
		{
		  double f = 1. - x/rad;
			return f < 0. ? 0. : f;
		}
	};
	//! Struct for a constant function, used in Filter2D and Filter3D
	struct constantFunction
	{
		const double val;
		explicit constantFunction(double inVal) : val(inVal) {}
		double operator()(double x) {return 1.;}
	};
	//! Struct to divide by a value
	struct divVal
	{
		double val;
		divVal(double inVal) : val(inVal){}
		double operator()(double inval){return inval/val;}
	};
	//! Struct to multiply by a value
	template<typename T>
	struct multVal
	{
		T val;
		multVal(const T& inVal) : val(inVal){}
		T operator()(const T& x) {return x*val;}
	};

	//! Struct to compare values, greater than
	struct greaterThanX
	{
		double x;
		greaterThanX(double inx) : x(inx) {}
		bool operator()(double val) {return (val > x);}
	};
	//! Struct to compare values, less than
	struct lessThanX
	{
		double x;
		lessThanX(double inx) : x(inx) {}
		bool operator()(double val) {return (val < x);}
	};
	//! Converts an int to a double
	inline double int2doub(int k){return (double)k;}
	//! Returns a bool indicating whether or not the argument is greater than 1
	inline bool greaterThan1 (double val) { return ( val > 1.); }
	//! Returns a bool indicating whether or not the argument is less than 0
	inline bool lessThan0 (double val) {return (val < 0.); }
	//! Norm of vector
	inline double norm(const std::vector<double>& v)
	{
		return sqrt(std::accumulate(v.begin(), v.end(), 0., HelperNS::powSum(2.)));
	}
	//! Sum two std::vector objects
	template<typename T> inline
	std::vector<T> vecSum(const std::vector<T>& v1, const std::vector<T>& v2)
	{
		std::vector<T> v3;
		assert(v1.size() == v2.size());
		v3.reserve(v1.size());
		std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(v3), std::plus<T>());
		return v3;
	}
	//! Subtract two std::vector objects
	template<typename T> inline
	std::vector<T> vecMinus(const std::vector<T>& v1, const std::vector<T>& v2)
	{
		std::vector<T> v3;
		assert(v1.size() == v2.size());
		v3.reserve(v1.size());
		std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(v3), std::minus<T>());
		return v3;
	}
	//! Scalar multiplication of a std::vector object
	template<typename T> inline
	std::vector<T> vecScalarMult(const std::vector<T>& v, const T& s)
	{
		std::vector<T> vout;
		vout.reserve(v.size());
		HelperNS::multVal<T> multOp(s);
		std::transform(v.begin(), v.end(), std::back_inserter(vout), multOp);
		return vout;
	}
	//! Scalar multiplication of a std::vector object
	template<typename T> inline
	std::vector<T> vecScalarMult(const T& s, const std::vector<T>& v)
	{
		return vecScalarMult(v, s);
	}
	//! Absolute difference between two std::size_t
	inline
	std::size_t absdiff(std::size_t st1, std::size_t st2)
	{
		return st1 > st2 ? st1 - st2 : st2 - st1;
	}
}
}
#endif

