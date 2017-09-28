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

#ifndef OUTPUTHANDLER_H
#define OUTPUTHANDLER_H

#include "inputloadertopopt.h"

namespace Topologies{
class TopOptRep;
class TopOptObjFun;

//! A class to handle writing output to a file
/*! This class handles output writing as specified by a set of input parameters.
 *  Output can be generated both at a specified number of optimization iterations and
 *  after completion.  There are various types of input including volume mesh, surface mesh,
 *  raw data, STL, etc.
 */
class OutputHandler
{
public:
	//! Constructor that takes a set of input paramters from the InputLoader
	OutputHandler(const InputLoader::Output& inputParams);
	~OutputHandler() {}

	//! Write output for `torToPrint` to file
	/*! This function should be called at each optimization iteration, as the number of calls
	 *  determines how often to write output.  Note then that calling this function will not
	 *  necessarily output anything, it depends on if it is time to output based on the specified
	 *  output period.  
	 */
	void handleOutput(const TopOptRep* const torToPrint, const TopOptObjFun* const toofFunc, bool lastOutput) const;

private:
	void dispatchOutput(const TopOptRep* const torToPrint, const TopOptObjFun* const toofFunc, const std::string& curFileName) const;
	void handleOutput2d(const TopOptRep* const torToPrint, const std::string& curFileName) const;
	void handleOutput3d(const TopOptRep* const torToPrint, const std::string& curFileName) const;
	void stripFileExtension();
	std::string getFileExtensionString() const;
private:
	mutable unsigned kout;
	mutable bool writtenOnce;
	bool overwrite;
	bool outputAtFin, outputStep;
	unsigned outputPeriod;
	OutputType type;
	std::string fileName;
	OutputFileFormat fileFormat;
	double extrusionLength;
};
}
#endif

