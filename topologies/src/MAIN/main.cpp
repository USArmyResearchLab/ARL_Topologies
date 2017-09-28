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

/*! @mainpage Topologies
 *
 * @section intro_sec Introduction
 *
 * Topologies is an extensible topology optimization program meant to be both a research
 * platform and hopefully usable in general.  To achieve this, it separates concepts from topology 
 * optimization into three parts: Representation, optimizer, and objective function.
 *
 * Topologies was originally developed at the US Army Research Laboratory, by Raymond Wildman.
 * Please contact raymond.a.wildman.civ@mail.mil with any questions.
 * 
 * @subsection repsec Representation
 *
 * A topology representation is the method in which a topology or geometry is parameterized
 * so that an optimization algorithm can optimize it.  The most common representation is a
 * pixel or voxel method wherein a region is broken into rectangular pixels or voxels and
 * an unknown quantity is defined on each.  An optimization algorithm then attempts to find
 * the optimal value for each pixel.  
 *
 * @sa TopOptRep
 *  
 * @subsection optsec Optimizer
 *
 * An optimization algorithm attempts to find a set of parameters that minimizes or
 * maximizes a given objective function.  There are a few main characteristics of 
 * optimization methods such as constrained/unconstrained and those that use gradients
 * and those that don't.  
 * 
 * @sa TopOpt
 *
 * @subsection objfunsec Objective Function
 * 
 * An objective function is the problem that the optimizer is attempting to minimize/maximize.
 * It is typically some type of physics problem implemented numerically, such as a finite element
 * solver.  No objective function is implemented directly into Topologies; objective functions
 * are implemented in a shared library that Topologies will load.  However, a basic finite 
 * element solver is included that solves the standard compliance minimization problem in
 * structural topology optimization.
 *
 * @sa TopOptObjFun
 *
 */

#include "topoptuniverse.h"

Topologies::VerbLevel verbosity = Topologies::vlSilent;

int main(int argc, char* argv[])
{
	Topologies::TopOptUniverse toProblem(argc, argv);
	toProblem.runProblem();
}

