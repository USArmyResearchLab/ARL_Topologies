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

#include "inputloaderopt.h"

namespace Topologies{
namespace InputLoader
{
	void TOOChain::parse(const pugi::xml_document& xmldoc)
	{
		parse(xmldoc.child("chain"));
	}

	void TOOChain::parse(const pugi::xml_node& rootNode)
	{
		OptNodeInfo tmpONI;
		for(pugi::xml_node stepNode = rootNode.child(tmpONI.getNodeName().c_str());
			stepNode; stepNode = stepNode.next_sibling(tmpONI.getNodeName().c_str()))
		{
			oniVec.push_back(tmpONI);
			oniVec.back().parse(stepNode, curFileName);
		}
	}

	void TOOGeneric::parse(const pugi::xml_document& xmldoc)
	{
		parse(xmldoc.child(optimizerName.c_str()));
	}

	void TOOGeneric::parse(const pugi::xml_node& rootNode)
	{
		// All optional (see defaults)
		try{maxStep = readDoublePCData(rootNode, "max_step");}
		catch(ParseException pe){}
		try{stepSize = readDoublePCData(rootNode, "step_size");}
		catch(ParseException pe){}
		try{filterSize = readDoublePCData(rootNode, "filter_size");}
		catch(ParseException pe){}
		try{constraintPenalty = readDoublePCData(rootNode, "constraint_penalty");}
		catch(ParseException pe){}
		try{penaltyPower = readDoublePCData(rootNode, "penalty_power");}
		catch(ParseException pe){}
		try{stopTol = readDoublePCData(rootNode, "stop_tol");}
		catch(ParseException pe){}
		try{maxIters = readUnsignedPCData(rootNode, "max_iterations");}
		catch(ParseException pe){}
	}
	
	void TOOGA::parse(const pugi::xml_document& xmldoc)
	{
		if(isPareto)
			parse(xmldoc.child("pareto_ga"));
		else
			parse(xmldoc.child("ga"));
	}

	void TOOGA::parse(const pugi::xml_node& rootNode)
	{
		// Required inputs
		try
		{
			numGens = readUnsignedPCData(rootNode, "num_generations");
			popSize = readUnsignedPCData(rootNode, "population_size");
			crossRate = readDoublePCData(rootNode, "crossover_rate");
			mutationRate = readDoublePCData(rootNode, "mutation_rate");
			mutationRange = readDoublePCData(rootNode, "mutation_range");
			if(isPareto)
			{
				numGoals = readUnsignedPCData(rootNode, "num_goals");
				goalWeights = readDoubleVecPCData(rootNode, "goal_weights");
				if(goalWeights.size() != numGoals)
				{
					std::string errorStr("goal_weights doesn't match num_goals");
					throw ParseException(petInvalidNumericalInput, std::move(errorStr));
				}
			}
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
		//Optional
		try{ntourn = readUnsignedPCData(rootNode, "num_tourn");}
		catch(ParseException pe){}
		try{numElite = readUnsignedPCData(rootNode, "num_elitism");}
		catch(ParseException pe){}
		try{mutationRadius = readUnsignedPCData(rootNode, "mutation_radius");}
		catch(ParseException pe){}
		try{constraintPenalty = readDoublePCData(rootNode, "constraint_penalty");}
		catch(ParseException pe){}
		try{penaltyPower = readDoublePCData(rootNode, "penalty_power");}
		catch(ParseException pe){}
		try{sharingRadius = readDoublePCData(rootNode, "sharing_radius");}
		catch(ParseException pe){}
	}
}
}

