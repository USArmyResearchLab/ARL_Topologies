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

 /*! @file topologiesdefs.h
  *	 @brief Enumerations, macros, and typedefs
 */

#ifndef TOPOLOGIESDEFS_H
#define TOPOLOGIESDEFS_H

#include <complex>

//#define USE_MPI // run in parallel with MPI
//#define USE_MULTICORE_GA // run chromosomes in parallel with MSVC concurrency
//#define USE_MULTICORE_LO // runs local optimizer gradient computations in parallel
//#define USE_MULTICORE_OBJFUN // run objective function in parallel with MSVC concurrency

#define MAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) 
#define TOL_EQ(a, b, tol) (fabs((a) - (b)) <= (tol))
#define TOL_LEQ(a, b, tol) ((a) < (b) || TOL_EQ(a, b, tol))
#define TOL_GEQ(a, b, tol) ((a) > (b) || TOL_EQ(a, b, tol))
#define TOL_LT(a, b, tol) ((a) < (b) && !TOL_EQ(a, b, tol))
#define TOL_GT(a, b, tol) ((a) > (b) && !TOL_EQ(a, b, tol))

namespace Topologies{
//! A convenience typdef to easily change the precision of the code
typedef double Real;
//! Complex type defined from Real
typedef std::complex<Real> Complex;
//! typedef of unsigned int
typedef unsigned int Uint;

//! The ratio of the circumference of a circle to its diameter.
const Real PI = 3.14159265358979323846264338327950;
//! A constant to convert degrees to radians
const Real Deg2Rad = PI/180.0L;
//! A constant to convert radians to degrees
const Real Rad2Deg = 1.L/Deg2Rad;

//! An enumeration to define level of text output to stdout
enum VerbLevel
{
	vlSilent = 0,
	vlLaconic,
	vlLoquacious,
	vlDebug
};

//! Enumeration for TopOptRep types
enum TORType{tortPixel, tortLowRankPixel, tortHeaviside2D, tortVoxel, tortCSG2D, tortCSG3D, tortAlpha2D, tortAlpha3D, tortMesh2D, tortMesh3D, tortHeavisideMesh2D, tortHeaviside3D, tortHeavisideMesh3D, tortLinearSpline, tortUnknown};
//! Enumeration for TopOpt optimizers
enum TOOType{tootOC, tootMMA, tootBFGS, tootGA, tootPGA, tootPSO, tootGD, tootChain, tootRefine, tootContinuation, tootConvert, tootUnknown};

//! An enumeration for defining temporal units
enum TimeUnits {tuTimeSteps, tuSeconds, tuMilliSeconds, tuMicroSeconds, tuNanoSeconds, tuPicoSeconds};
//! An enumeration for defining spatial units
enum SpaceUnits {suXSteps, suYSteps, suZSteps, suKilometers, suMeters, suCentimeters, suMillimeters, suMicrometers, suNanometers, suPicometers};
//! Function evaluation type, either objective function, constraint, gradient, or constraint gradient
enum EvalFunction{efF, efC, efG, efGC};
//! A data structure holding discretization parameters for space and time, and whether or not they have been set
struct DiscretizationParameters
{
	DiscretizationParameters() : dxSet(false), dySet(false), dzSet(false), dtSet(false), dx(1), dy(1), dz(1), dt(1){}
	bool dxSet, dySet, dzSet, dtSet;
	Real dx, dy, dz, dt;
};
}
#endif

