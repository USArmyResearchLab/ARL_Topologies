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

#include "outputhandler.h"
#include "topoptrep.h"
#include "outputwriter.h"
#include "postprocess.h"
#include "topoptobjfun.h"
#include "tomesh.h"
#include "mpihandler.h"
#include <memory>

namespace Topologies{
OutputHandler::OutputHandler(const InputLoader::Output& inputParams) :
	type(inputParams.getOutputType()),
	fileName(inputParams.getFileName()),
	fileFormat(inputParams.getFileFormat()),
	extrusionLength(inputParams.getExtrusionLength()),
	overwrite(inputParams.getOverwrite()),
	outputAtFin(inputParams.getOutputFinal()),
	outputStep(inputParams.getOutputPeriodic()),
	outputPeriod(inputParams.getOutputPeriod()),
	kout(0),
	writtenOnce(false)
{
	stripFileExtension();
	if(fileFormat == offGMSH)
		overwrite = true; // Add all output to 1 file
}

void OutputHandler::handleOutput(const TopOptRep* const torToPrint, const TopOptObjFun* const toofFunc, bool lastOutput, 
	MPIHandler* const mpih) const
{
	// Generate appropriate file name
	if(lastOutput && outputAtFin)
	{
		std::string curFileName = fileName;
		if(!overwrite)
			curFileName += "_fin" + getFileExtensionString();
		else
			curFileName += getFileExtensionString();
		// Generate output file
		dispatchOutput(torToPrint, toofFunc, mpih, curFileName);
	}
	else if(!lastOutput && outputStep && (++kout % outputPeriod) == 0)
	{
		std::string curFileName = fileName;
		if(!overwrite)
			curFileName += std::to_string(kout/outputPeriod) + getFileExtensionString();
		else
			curFileName += getFileExtensionString();
		// Generate output file
		dispatchOutput(torToPrint, toofFunc, mpih, curFileName);
	}
}

void OutputHandler::dispatchOutput(const TopOptRep* const torToPrint, const TopOptObjFun* const toofFunc, 
																	MPIHandler* const mpih, const std::string& curFileName) const
{
	if(type == otObjFunRes)
	{
		if(mpih == nullptr)
			toofFunc->printResult(*torToPrint, curFileName);
		else
			mpih->rootEvaluateTORAndPrint(torToPrint, curFileName);
	}
	else
	{
		if(torToPrint->getDimension() == 2)
			handleOutput2d(torToPrint, curFileName);
		else
			handleOutput3d(torToPrint, curFileName);
	}
	writtenOnce = true;
}

void OutputHandler::handleOutput2d(const TopOptRep* const torToPrint, const std::string& curFileName) const
{
	if(type == otSurface)
	{
		std::vector<Mesh_Segment_2> segVec;
		torToPrint->get2DSegments(segVec);
		OutputWriter::plotSegmentsMatlab(segVec, curFileName);
		if(fileFormat != offMatlab && fileFormat != offDefault)
			std::cout << "Warning: can't output 2D segments to STL or VTK, using Matlab for file: " << fileName << std::endl;
	}
	else if(type == otVolume)
	{
		std::unique_ptr<TOMesh> outMesh(torToPrint->getOutputMesh());
		if(fileFormat == offMatlab)
			OutputWriter::plotMeshMatlab(outMesh.get(), curFileName);
		else if(fileFormat == offVTK)
			OutputWriter::plotMeshVTK(outMesh.get(), curFileName);
		else if(fileFormat == offGMSH)
			OutputWriter::plotMeshGMSH(outMesh.get(), curFileName, kout, !writtenOnce);
		else if(fileFormat == offSTL)
			std::cout << "Warning: can't output 2D mesh to STL, using Matlab for file: " << fileName << std::endl;
	}
	else if(type == otExtrude)
	{
		PostProcess::linearExtrude(*torToPrint, curFileName, extrusionLength);
		if(fileFormat != offSTL)
			std::cout << "Warning: Extrusion output only supported using STL output, using STL for file: " << fileName << std::endl;
	}
	else if(type == otRawData)
	{
		if(fileFormat == offMatlab)
		{
			std::vector<double> realVec;
			torToPrint->getRealRep(realVec);
			std::vector<std::size_t> dataSizes;
			torToPrint->getDataSize(dataSizes);
			if(dataSizes.size() >= 3)
				OutputWriter::plotDataMatlab(realVec, dataSizes[0], dataSizes[1], dataSizes[2], curFileName);
			else
				std::cout << "Error writing raw_data file " << fileName << ", didn't find correct data sizes" << std::endl;
		}
		else
			OutputWriter::outputTORData(torToPrint, curFileName);
	}
}

void OutputHandler::handleOutput3d(const TopOptRep* const torToPrint, const std::string& curFileName) const
{
	if(type == otSurface)
	{
		std::unique_ptr<TOMesh> outMesh(torToPrint->get3DSurfaceMesh());
		OutputWriter::plotSTL(outMesh.get(), curFileName);
		if(fileFormat == offMatlab)
			std::cout << "Warning: 3d surface output to Matlab not currently supported, using STL for file: " << fileName << std::endl;
	}
	else if(type == otVolume)
	{
		std::unique_ptr<TOMesh> outMesh(torToPrint->getOutputMesh());
		if(fileFormat == offMatlab)
			OutputWriter::plotMeshMatlab(outMesh.get(), curFileName);
		else if(fileFormat == offVTK)
			OutputWriter::plotMeshVTK(outMesh.get(), curFileName);
		else if(fileFormat == offGMSH)
			OutputWriter::plotMeshGMSH(outMesh.get(), curFileName, kout, !writtenOnce);
		else if(fileFormat == offSTL)
			std::cout << "Warning: can't output 3D mesh to STL, using Matlab for file: " << fileName << std::endl;
	}
	else if(type == otExtrude)
		std::cout << "Warning in Output: Can't extrude a 3d design, no output for file: " << fileName << std::endl;
	else if(type == otRawData)
	{
		if(fileFormat == offMatlab)
		{
			std::vector<double> realVec;
			torToPrint->getRealRep(realVec);
			std::vector<std::size_t> dataSizes;
			torToPrint->getDataSize(dataSizes);
			if(dataSizes.size() >= 3)
				OutputWriter::plotDataMatlab(realVec, dataSizes[0], dataSizes[1], dataSizes[2], curFileName);
			else
				std::cout << "Error writing raw_data file " << fileName << ", didn't find correct data sizes" << std::endl;
		}
		else
			OutputWriter::outputTORData(torToPrint, curFileName);
	}
}

std::string OutputHandler::getFileExtensionString() const
{
	if(fileFormat == offMatlab)
		return ".m";
	else if(fileFormat == offSTL)
		return ".stl";
	else if(fileFormat == offVTK)
		return ".vtk";
	else if(fileFormat == offGMSH)
		return ".msh";
	return "";
}

void OutputHandler::stripFileExtension()
{
	std::string extension = getFileExtensionString();
	if(!extension.empty())
	{
		std::size_t n = fileName.rfind(extension);
		if(n != std::string::npos) // Extension exists in file name
			fileName.replace(n, extension.length(), ""); // Remove extension
	}
}
}

