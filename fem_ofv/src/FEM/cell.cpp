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

#include "cell.h"

Cell::Cell(const Topologies::GenericMaterial& inMat) :
	itsMaterial(inMat)
{
}

unsigned short int Cell::getLocalBF(unsigned BFNum) const
{
	for (unsigned short int i = 0; i < bfVec.size(); ++i)
		if(BFNum == bfVec[i])
			return i;
	return bfVec.size();
}

