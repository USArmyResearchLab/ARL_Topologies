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

#include "csgnode.h"

namespace Topologies{
using std::vector;

vector<GenericMaterial> CSGNode::materialList;
DecodeType CSGNode::theDT = dtBoolean;
DimensionType CSGNode::theDim = dt2d;
int CSGNode::numFuncTypes = 2;
Real CSGNode::boundMag = 1, CSGNode::boundSnap = 0.1, CSGNode::prioritymax = 10.;;
bool CSGNode::useBoundSnap = false;

CSGNode::CSGNode(int curTreeDepth):
	parent(nullptr),
	depthInTree(curTreeDepth + 1)
{
}

CSGNode::CSGNode(const CSGNode& copy) :
	parent(nullptr),
	depthInTree(copy.depthInTree)
{
}

CSGNode::CSGNode(CSGNode && copy)
{
	swap(copy);
}

void CSGNode::swap(CSGNode& arg2)
{
	std::swap(depthInTree, arg2.depthInTree);
}

CSGNode::~CSGNode()
{
}
}

