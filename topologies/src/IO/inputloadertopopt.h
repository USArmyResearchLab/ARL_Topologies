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

#ifndef INPUTLOADERTOPOPT_H
#define INPUTLOADERTOPOPT_H

#include "inputloader.h"

namespace Topologies{
class TopOptRep;

namespace InputLoader
{
	//! A class for parsing the input file for output settings
	class Output : public InputParser
	{
	public:
		Output() {}
		virtual ~Output() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns the OutputType
		OutputType getOutputType() const {return type;}
		//! Returns the file name to write output
		const std::string& getFileName() const {return fileName;}
		//! Returns the OutputFileFormat
		OutputFileFormat getFileFormat() const {return fileFormat;}
		//! Returns the extrusion length
		double getExtrusionLength() const {return extrusionLength;}
		//! Returns the output period
		unsigned getOutputPeriod() const {return outputPeriod;}
		//! Returns whether or not to overwrite output files
		bool getOverwrite() const {return overwrite;}
		//! Returns whether or not to write output after completion
		bool getOutputFinal() const {return outputAtFin;}
		//! Returns whether or not to write output periodically (according to outputPeriod)
		bool getOutputPeriodic() const {return outputStep;}
	private:
		void setOutputType(const pugi::xml_node& rootNode);
		void setOutputFileFormat(const pugi::xml_node& rootNode);

		OutputType type;
		std::string fileName;
		OutputFileFormat fileFormat;
		double extrusionLength;
		unsigned outputPeriod = 1;
		bool overwrite = true, outputAtFin = true, outputStep = false;
	};

	//! A class for loading the initial guess to use for the optimizer
	class InitialGuess : public InputParser
	{
	public:
		InitialGuess() {}
		virtual ~InitialGuess() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns the InitialGuessType
		InitialGuessType getInitialGuessType() const {return theIGT;}
		//! Returns the file name of the initial guess to load
		const std::string& getFileName() const {return fileName;}
		//! Returns the initial guess constant
		double getConstant() const {return igConstant;}
		//! Returns the initial guess range for random number generation (uniform distribution)
		const std::pair<double,double> getRandRange() const {return randRange;}
	private:
		void setIGT(const pugi::xml_node& rootNode);
		bool setIGTConst(const pugi::xml_node& rootNode);
		bool setIGTRandMin(const pugi::xml_node& rootNode);
		bool setIGTRandMax(const pugi::xml_node& rootNode);

		InitialGuessType theIGT = igtConstant;
		std::string fileName;
		double igConstant = 0.5;
		std::pair<double, double> randRange = std::pair<double,double>(0., 1.);
	};

	//! A class for loading the main functionality of the code: Representation, optimizer, and objective function
	/*! This class handles the input at the top level.  This will parse the main input file to get the basic
	*  information needed to continue input parsing.  It loads the TORep type, TOOptimizer type, and the objective function
	*  input file name and shared library name.  Parsing is then passed to the files responsible for those other inputs.
	*/
	class TopologiesParser : public InputParser
	{
	public:
		TopologiesParser() {}
		virtual ~TopologiesParser() {}
		virtual void parse(const pugi::xml_document& xmldoc);
		virtual void parse(const pugi::xml_node& rootNode);

		//! Returns the file name of the objective function input
		std::string getTOFFileName() const {return tofInputFileName;}
		//! Returns the file name of the objective function shared library
		std::string getTOFSharedLibName() const {return tofSharedLibFileName;}
		//! Returns the number of processors to use per objective function evaluation
		unsigned getNumProcessorsPerOFV() const {return numProcessorsPerOFV;}
		//! Returns the InitialGuess object that contains input values
		const InitialGuess& getInitialGuess() const {return igParser;}
		//! Returns an iterator to the first Output object
		std::vector<Output>::const_iterator outputBegin() const {return opVec.begin();}
		//! Returns the past-the-end iterator of the Output object vector
		std::vector<Output>::const_iterator outputEnd() const {return opVec.end();}
		//! Returns the RepNodeInfo object holding the location of the TopOptRep to load
		const RepNodeInfo& getRepNodeInfo() const {return mainRNI;}
		//! Returns the OptNodeInfo object holding the location of the TopOpt to load
		const OptNodeInfo& getOptNodeInfo() const {return mainONI;}
	private:
		std::string tofInputFileName, tofSharedLibFileName;
		RepNodeInfo mainRNI;
		OptNodeInfo mainONI;
		unsigned numProcessorsPerOFV = 1;

		InitialGuess igParser;
		std::vector<Output> opVec;
	};
}
}
#endif

