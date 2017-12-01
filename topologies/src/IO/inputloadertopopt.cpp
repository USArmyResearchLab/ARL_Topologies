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

#include "inputloadertopopt.h"
#include "helper.h"
#include "topoptrep.h"
#include <string>
#include <fstream>

namespace Topologies{
namespace InputLoader
{
	void TopologiesParser::parse(const pugi::xml_document& xmldoc)
	{
		parse(xmldoc.child("topologies"));
	}

	void TopologiesParser::parse(const pugi::xml_node& rootNode)
	{
		// First parse required inputs
		mainRNI.parse(rootNode.child(mainRNI.getNodeName().c_str()), curFileName);
		mainONI.parse(rootNode.child(mainONI.getNodeName().c_str()), curFileName);
		try
		{
			// Parse objective function
			tofInputFileName = readAndCheckFileNameAttribute(rootNode, "objective_function", "input_file");
			tofSharedLibFileName = readStringAttribute(rootNode, "objective_function", "shared_library");
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
		// Now parse unrequired options
		try{numProcessorsPerOFV = readUnsignedPCData(rootNode, "num_processors_per_ofv");}
		catch(ParseException pe){} // Not required

		// Parse initial guess
		pugi::xml_node igNode = rootNode.child("initial_guess");
		if(igNode) // InitialGuess defaults to a constant value of 0.5
			igParser.parse(igNode);

		// Parser output
		for(pugi::xml_node cn = rootNode.child("output"); cn; cn = cn.next_sibling("output"))
		{
			opVec.push_back(Output());
			opVec.back().parse(cn);
		}
	}

	void InitialGuess::parse(const pugi::xml_document& xmldoc)
	{
		parse(xmldoc.child("topologies").child("initial_guess"));
	}

	void InitialGuess::parse(const pugi::xml_node& rootNode)
	{
		// Read required inputs
		try
		{
			setIGT(rootNode);
			// Conditionally required inputs
			if(theIGT == igtFile)
				fileName = readAndCheckFileNamePCData(rootNode, "file_name");
			if(theIGT == igtConstant || theIGT == igtConstantWithNoise)
				igConstant = readDoublePCData(rootNode, "constant_val");
			if(theIGT == igtRandom || theIGT == igtConstantWithNoise)
			{
				randRange.first = readDoublePCData(rootNode, "noise_min");
				randRange.second = readDoublePCData(rootNode, "noise_max");
			}
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
	}

	void InitialGuess::setIGT(const pugi::xml_node& rootNode)
	{
		std::string retval = readStringAttribute(rootNode, "type");
		theIGT = parserIGType(retval);
		if(theIGT == igtUnknown)
			throw ParseException(petUnknownInput, std::move(retval));
	}

	void Output::parse(const pugi::xml_document& xmldoc)
	{
		 parse(xmldoc.child("topologies").child("output"));
	}

	void Output::parse(const pugi::xml_node& rootNode)
	{
		// Read required inputs
		try
		{
			setOutputType(rootNode);
			setOutputFileFormat(rootNode);
			fileName = readStringPCData(rootNode, "file_name");
			// Conditionally optional
			if(type == otExtrude)
				extrusionLength = readDoublePCData(rootNode, "extrusion_length");
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
		// Optional inputs
		try{overwrite = readBoolPCData(rootNode, "overwrite");}
		catch(ParseException pe){}
		try{outputPeriod = readUnsignedPCData(rootNode, "output_period");}
		catch(ParseException pe){}
		try{outputAtFin = readBoolPCData(rootNode, "write_final_result");}
		catch(ParseException pe){}
		try{outputStep = readBoolPCData(rootNode, "write_periodic_results");}
		catch(ParseException pe){}
	}

	void Output::setOutputType(const pugi::xml_node& rootNode)
	{
		std::string retval = readStringAttribute(rootNode, "type");
		type = parseOutputType(retval);
		if(type == otUnknown)
			throw ParseException(petUnknownInput, std::move(retval));
	}

	void Output::setOutputFileFormat(const pugi::xml_node& rootNode)
	{
		std::string retval = readStringPCData(rootNode, "file_format");
		fileFormat = parseOutputFileFormat(retval);
	}
}
}

