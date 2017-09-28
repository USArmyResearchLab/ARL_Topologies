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

#include "inputloader.h"
#include "helper.h"
#include "topoptrep.h"
#include <string>
#include <fstream>

namespace Topologies{
VolMeshTORSpecification::VolMeshTORSpecification() :
	torDimType(dt2d), torMeshType(mffStructured), torUnknownLocation(ulElement), torFilterType(ftNone),
	torPenalizationType(ptSIMP), torProjectionType(ptNone)
{
}

VolMeshTORSpecification::VolMeshTORSpecification(DimensionType inDT, MeshFileFormat inMFF, UnknownLocation inUL, FilterType inFT,
	PenalizationType inPenType, ProjectionType inProjType) :
		torDimType(inDT), torMeshType(inMFF), torUnknownLocation(inUL), torFilterType(inFT),
		torPenalizationType(inPenType), torProjectionType(inProjType)
{
}

VolMeshTORSpecification::VolMeshTORSpecification(TORType inTORT)
{
	if(inTORT == tortPixel || inTORT == tortVoxel || inTORT == tortMesh2D || inTORT == tortMesh3D)
	{
		torDimType = inTORT == tortPixel || inTORT == tortMesh2D? dt2d : dt3d;
		torMeshType = inTORT == tortPixel || inTORT == tortVoxel ? mffStructured : mffUnknown;
		torUnknownLocation = ulElement;
		torFilterType = ftNone;
		torPenalizationType = ptSIMP;
		torProjectionType = ptNone;
	}
	else if(inTORT == tortHeaviside2D || inTORT == tortHeaviside3D || inTORT == tortHeavisideMesh2D || inTORT == tortHeavisideMesh3D)
  {
    torDimType = inTORT == tortHeaviside2D || inTORT == tortHeavisideMesh2D ? dt2d : dt3d;
    torMeshType = inTORT == tortHeaviside2D || inTORT == tortHeaviside3D ? mffStructured : mffUnknown;
    torUnknownLocation = ulNode;
    torFilterType = ftLinear;
    torPenalizationType = ptSIMP;
    torProjectionType = ptThresholdHeavi;
  }
}

VolMeshTORSpecification::VolMeshTORSpecification(const std::vector<int>& params) : 
	torDimType((DimensionType)params[0]), torMeshType((MeshFileFormat)params[1]), torUnknownLocation((UnknownLocation)params[2]), 
	torFilterType((FilterType)params[3]), torPenalizationType((PenalizationType)params[4]), 
	torProjectionType((ProjectionType)params[5])
{
	assert(params.size() == 6);
}

std::vector<int> VolMeshTORSpecification::toVec() const 
{
	return {torDimType, torMeshType, torUnknownLocation, torFilterType, torPenalizationType, torProjectionType};
}

std::size_t VolMeshTORSpecification::size() const
{
	return toVec().size();
}

namespace InputLoader
{
	void inputTORData(TopOptRep* const inTOR, const std::string& fileName)
	{
		std::ifstream inFile(fileName);
		std::vector<std::vector<int>> discreteVars;
		std::size_t curSize;
		inFile >> curSize;
		discreteVars.resize(curSize);
		for(std::size_t k1 = 0; k1 < discreteVars.size(); ++k1)
		{
			inFile >> curSize;
			discreteVars[k1].resize(curSize);
			for(std::size_t k2 = 0; k2 < discreteVars[k1].size(); ++k2)
				inFile >> discreteVars[k1][k2];
		}
		std::vector<std::vector<double>> realVars;
		inFile >> curSize;
		realVars.resize(curSize);
		for(std::size_t k1 = 0; k1 < realVars.size(); ++k1)
		{
			inFile >> curSize;
			realVars[k1].resize(curSize);
			for(std::size_t k2 = 0; k2 < realVars[k1].size(); ++k2)
				inFile >> realVars[k1][k2];
		}
		inTOR->setMPIRep(discreteVars, realVars);
	}

	InitialGuessType parserIGType(const std::string& strIGType)
	{
		if(HelperNS::upperCase(strIGType) == "RANDOM")
			return igtRandom;
		else if(HelperNS::upperCase(strIGType) == "CONSTANT")
			return igtConstant;
		else if(HelperNS::upperCase(strIGType) == "CONSTANT_WITH_NOISE")
			return igtConstantWithNoise;
		else if(HelperNS::upperCase(strIGType) == "FILE")
			return igtFile;
		return igtUnknown;
	}

	TORType parseTORType(const std::string& strTORType)
	{
		//tortPixel, tortLowRankPixel, tortVoxel, tortCSG2D, tortCSG3D, tortAlpha2D, tortAlpha3D
		if(HelperNS::upperCase(strTORType) == "PIXEL")
			return tortPixel;
		else if(HelperNS::upperCase(strTORType) == "HEAVISIDE2D")
			return tortHeaviside2D;
		else if(HelperNS::upperCase(strTORType) == "LOW_RANK_PIXEL")
			return tortLowRankPixel;
		else if(HelperNS::upperCase(strTORType) == "VOXEL")
			return tortVoxel;
		else if(HelperNS::upperCase(strTORType) == "HEAVISIDE3D")
			return tortHeaviside3D;
		else if(HelperNS::upperCase(strTORType) == "CSG2D")
			return tortCSG2D;
		else if(HelperNS::upperCase(strTORType) == "CSG3D")
			return tortCSG3D;
		else if(HelperNS::upperCase(strTORType) == "ALPHA2D")
			return tortAlpha2D;
		else if(HelperNS::upperCase(strTORType) == "ALPHA3D")
			return tortAlpha3D;
		else if(HelperNS::upperCase(strTORType) == "MESH2D")
			return tortMesh2D;
		else if(HelperNS::upperCase(strTORType) == "MESH3D")
			return tortMesh3D;
		else if(HelperNS::upperCase(strTORType) == "HEAVISIDEMESH2D")
      return tortHeavisideMesh2D;
    else if(HelperNS::upperCase(strTORType) == "HEAVISIDEMESH3D")
      return tortHeavisideMesh3D;
		return tortUnknown;
	}

	TOOType parseTOOType(const std::string& strTOOType)
	{
		//tootOC, tootLocal, tootGA, tootPSO, tootUnknown
		if(HelperNS::upperCase(strTOOType) == "OC")
			return tootOC;
		else if(HelperNS::upperCase(strTOOType) == "MMA")
			return tootMMA;
		else if(HelperNS::upperCase(strTOOType) == "BFGS")
			return tootBFGS;
		else if(HelperNS::upperCase(strTOOType) == "GA")
			return tootGA;
		else if(HelperNS::upperCase(strTOOType) == "PGA")
			return tootPGA;
		else if(HelperNS::upperCase(strTOOType) == "PSO")
			return tootPSO;
		else if(HelperNS::upperCase(strTOOType) == "GD")
			return tootGD;
		else if(HelperNS::upperCase(strTOOType) == "CHAIN")
			return tootChain;
		else if(HelperNS::upperCase(strTOOType) == "REFINE")
			return tootRefine;
		else if(HelperNS::upperCase(strTOOType) == "CONVERT")
			return tootConvert;
		else if(HelperNS::upperCase(strTOOType) == "CONTINUATION")
			return tootContinuation;
		return tootUnknown;
	}

	OutputType parseOutputType(const std::string& strOType)
	{
		if(HelperNS::upperCase(strOType) == "SURFACE")
			return otSurface;
		else if(HelperNS::upperCase(strOType) == "VOLUME")
			return otVolume;
		else if(HelperNS::upperCase(strOType) == "EXTRUDE")
			return otExtrude;
		else if(HelperNS::upperCase(strOType) == "OBJECTIVE_FUNCTION_RESULT")
			return otObjFunRes;
		else if(HelperNS::upperCase(strOType) == "RAW_DATA")
			return otRawData;
		return otUnknown;
	}

	OutputFileFormat parseOutputFileFormat(const std::string& strOFF)
	{
		if(HelperNS::upperCase(strOFF) == "MATLAB")
			return offMatlab;
		else if(HelperNS::upperCase(strOFF) == "STL")
			return offSTL;
		else if(HelperNS::upperCase(strOFF) == "VTK")
			return offVTK;
		else if(HelperNS::upperCase(strOFF) == "GMSH")
			return offGMSH;
		else if(HelperNS::upperCase(strOFF) == "DEFAULT")
			return offDefault;
		return offDefault;
	}

	MeshElementType parseMeshElementType(const std::string& strMEType)
	{
		if(HelperNS::upperCase(strMEType) == "TRI")
			return metTri;
		else if(HelperNS::upperCase(strMEType) == "QUAD")
			return metQuad;
		else if(HelperNS::upperCase(strMEType) == "TET")
			return metTet;
		else if(HelperNS::upperCase(strMEType) == "HEX")
			return metHex;
		return metUnknown;
	}

	MeshFileFormat parseMeshFileFormat(const std::string& strMFF)
	{
		if(HelperNS::upperCase(strMFF) == "EXODUSII")
			return mffExodus;
		else if(HelperNS::upperCase(strMFF) == "GMSH")
			return mffGMSH;
		else if(HelperNS::upperCase(strMFF) == "POLYGON")
			return mffPolygon;
		else if(HelperNS::upperCase(strMFF) == "STL")
			return mffSTL;
		return mffUnknown;
	}

	namespace
	{
		pugi::xml_node traversePCDataNode(const pugi::xml_node& rootNode, const std::vector<std::string>& path)
		{
			std::string fullPath;
			pugi::xml_node targetNode = rootNode;
			for(std::size_t k = 0; k < path.size(); ++k)
			{
				fullPath += path[k];
				fullPath += ":";
				targetNode = targetNode.child(path[k].c_str());
			}
			targetNode = targetNode.first_child();
			if(!targetNode)
				throw ParseException(petPCDataNotFound, std::move(fullPath));
			else if(targetNode.type() != pugi::node_pcdata)
				throw ParseException(petWrongNodeType, std::string(fullPath));
			return targetNode;
			
		}
		pugi::xml_attribute traverseAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName)
		{
			std::string fullPath;
			pugi::xml_node targetNode = rootNode;
			for(std::size_t k = 0; k < path.size(); ++k)
			{
				fullPath += path[k];
				fullPath += ":";
				targetNode = targetNode.child(path[k].c_str());
			}
			fullPath += attrName;
			pugi::xml_attribute targetAttr = targetNode.attribute(attrName.c_str());
			if(!targetAttr)
				throw ParseException(petAttributeNotFound, std::move(fullPath));
			return targetAttr;
		}
	}

	std::string readStringPCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path)
	{
		pugi::xml_node pcdataNode = traversePCDataNode(rootNode, path);
		return std::string(pcdataNode.value());
	}

	std::string readStringAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName)
	{
		pugi::xml_attribute attr = traverseAttribute(rootNode, path, attrName);
		return std::string(attr.value());
	}

	double readDoublePCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path)
	{
		pugi::xml_node pcdataNode = traversePCDataNode(rootNode, path);
		return pcdataNode.text().as_double();
	}

	double readDoubleAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName)
	{
		pugi::xml_attribute attr = traverseAttribute(rootNode, path, attrName);
		return attr.as_double();
	}

	unsigned readUnsignedPCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path)
	{
		pugi::xml_node pcdataNode = traversePCDataNode(rootNode, path);
		return pcdataNode.text().as_uint();
	}

	unsigned readUnsignedAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName)
	{
		pugi::xml_attribute attr = traverseAttribute(rootNode, path, attrName);
		return attr.as_uint();
	}

	int readIntPCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path)
	{
		pugi::xml_node pcdataNode = traversePCDataNode(rootNode, path);
		return pcdataNode.text().as_int();
	}

	int readIntAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName)
	{
		pugi::xml_attribute attr = traverseAttribute(rootNode, path, attrName);
		return attr.as_int();
	}

	bool readBoolPCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path)
	{
		pugi::xml_node pcdataNode = traversePCDataNode(rootNode, path);
		return pcdataNode.text().as_bool();
	}

	bool readBoolAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName)
	{
		pugi::xml_attribute attr = traverseAttribute(rootNode, path, attrName);
		return attr.as_bool();
	}

	std::string readStringPCData(const pugi::xml_node& rootNode, const std::string& path)
	{
		return readStringPCData(rootNode, std::vector<std::string>(1,path));
	}

	std::string readStringAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName)
	{
		return readStringAttribute(rootNode, std::vector<std::string>(1,path), attrName);
	}

	double readDoublePCData(const pugi::xml_node& rootNode, const std::string& path)
	{
		return readDoublePCData(rootNode, std::vector<std::string>(1,path));
	}

	double readDoubleAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName)
	{
		return readDoubleAttribute(rootNode, std::vector<std::string>(1,path), attrName);	
	}

	unsigned readUnsignedPCData(const pugi::xml_node& rootNode, const std::string& path)
	{
		return readUnsignedPCData(rootNode, std::vector<std::string>(1,path));	
	}

	unsigned readUnsignedAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName)
	{
		return readUnsignedAttribute(rootNode, std::vector<std::string>(1,path), attrName);
	}

	int readIntPCData(const pugi::xml_node& rootNode, const std::string& path)
	{
		return readIntPCData(rootNode, std::vector<std::string>(1,path));
	}

	int readIntAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName)
	{
		return readIntAttribute(rootNode, std::vector<std::string>(1,path), attrName);
	}

	bool readBoolPCData(const pugi::xml_node& rootNode, const std::string& path)
	{
		return readBoolPCData(rootNode, std::vector<std::string>(1,path));
	}

	bool readBoolAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName)
	{
		return readBoolAttribute(rootNode, std::vector<std::string>(1,path), attrName);
	}

	std::string readStringAttribute(const pugi::xml_node& rootNode, const std::string& attrName)
	{
		return readStringAttribute(rootNode, std::vector<std::string>(), attrName);
	}

	double readDoubleAttribute(const pugi::xml_node& rootNode, const std::string& attrName)
	{
		return readDoubleAttribute(rootNode, std::vector<std::string>(), attrName);	
	}

	unsigned readUnsignedAttribute(const pugi::xml_node& rootNode, const std::string& attrName)
	{
		return readUnsignedAttribute(rootNode, std::vector<std::string>(), attrName);
	}

	int readIntAttribute(const pugi::xml_node& rootNode, const std::string& attrName)
	{
		return readIntAttribute(rootNode, std::vector<std::string>(), attrName);
	}

	bool readBoolAttribute(const pugi::xml_node& rootNode, const std::string& attrName)
	{
		return readBoolAttribute(rootNode, std::vector<std::string>(), attrName);
	}

	std::string readStringPCData(const pugi::xml_node& rootNode)
	{
		return readStringPCData(rootNode, std::vector<std::string>());
	}

	double readDoublePCData(const pugi::xml_node& rootNode)
	{
		return readDoublePCData(rootNode, std::vector<std::string>());
	}

	unsigned readUnsignedPCData(const pugi::xml_node& rootNode)
	{
		return readUnsignedPCData(rootNode, std::vector<std::string>());
	}

	int readIntPCData(const pugi::xml_node& rootNode)
	{
		return readIntPCData(rootNode, std::vector<std::string>());
	}

	bool readBoolPCData(const pugi::xml_node& rootNode)
	{
		return readBoolPCData(rootNode, std::vector<std::string>());
	}

	std::vector<double> readDoubleVecPCData(const pugi::xml_node& rootNode, const std::string& childName)
	{
		std::vector<double> retval;
		std::string locname = childName;
		pugi::xml_node pcdNode = rootNode.child(childName.c_str());
		if(!pcdNode)
			throw ParseException(petPCDataNotFound, std::move(locname));
		for(pugi::xml_node stepNode = pcdNode.first_child();
					stepNode; stepNode = stepNode.next_sibling())
		{
			if(stepNode.type() == pugi::node_pcdata)
				retval.push_back(stepNode.text().as_double());
		}
		return retval;
	}
	std::vector<unsigned> readUnsignedVecPCData(const pugi::xml_node& rootNode, const std::string& childName)
	{
		std::vector<unsigned> retval;
		std::string locname = childName;
		pugi::xml_node pcdNode = rootNode.child(childName.c_str());
		if(!pcdNode)
			throw ParseException(petPCDataNotFound, std::move(locname));
		for(pugi::xml_node stepNode = pcdNode.first_child();
					stepNode; stepNode = stepNode.next_sibling())
		{
			if(stepNode.type() == pugi::node_pcdata)
				retval.push_back(stepNode.text().as_uint());
		}
		return retval;
	}
	std::vector<int> readIntVecPCData(const pugi::xml_node& rootNode, const std::string& childName)
	{
		std::vector<int> retval;
		std::string locname = childName;
		pugi::xml_node pcdNode = rootNode.child(childName.c_str());
		if(!pcdNode)
			throw ParseException(petPCDataNotFound, std::move(locname));
		for(pugi::xml_node stepNode = pcdNode.first_child();
				stepNode; stepNode = stepNode.next_sibling())
		{
			if(stepNode.type() == pugi::node_pcdata)
				retval.push_back(stepNode.text().as_int());
		}
		return retval;
	}
	std::vector<std::string> getNodePath(const pugi::xml_node& rootNode)
	{
		std::vector<std::string> path;
		for(pugi::xml_node cn = rootNode; cn; cn = cn.parent())
			path.push_back(std::string(cn.name()));
		path.pop_back(); // Top level node (it's nameless)
		std::reverse(path.begin(), path.end());
		return path;
	}

	std::string getExceptionName(const ParseExceptionType& inPE)
	{
		if(inPE == petPCDataNotFound)
			return "XML PCData node";
		if(inPE == petAttributeNotFound)
			return "XML Attribute";
		if(inPE == petUnknownInput)
			return "Unknown input";
		if(inPE == petInvalidNumericalInput)
			return "Invalid numerical input";
		if(inPE == petWrongNodeType)
			return "Incorrect node type encountered";
		return "Unknown exception!";
	}

	void errorMessage(const ParseException& inPE, bool willExit)
	{
		std::cerr << "Error in parsing input file : " << std::endl;
		std::cerr << "  " << getExceptionName(inPE.pet) << std::endl;
		std::cerr << "  " << inPE.info << std::endl;
		if(willExit)
			abort();
	}

	void NodeInfo::parse(const pugi::xml_node& rootNode, const std::string& curFileName)
	{
		// Required inputs
		try
		{
			tag = readStringAttribute(rootNode, "tag");
			path = getNodePath(rootNode);
			parseType(rootNode);
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
		// File name is optional, defaults to current file
		try{fileName = readStringAttribute(rootNode, "input_file");}
		catch(ParseException pe){fileName = curFileName;} // Set to current file
	}

	void RepNodeInfo::parseType(const pugi::xml_node& rootNode)
	{
		std::string retval = readStringAttribute(rootNode, "type");
		name = retval;
		nodeTORT = parseTORType(retval);
		if(nodeTORT == tortUnknown)
			throw ParseException(petUnknownInput, std::move(retval));
	}

	void OptNodeInfo::parseType(const pugi::xml_node& rootNode)
	{
		std::string retval = readStringAttribute(rootNode, "type");
		name = retval;
		nodeTOOT = parseTOOType(retval);
		if(nodeTOOT == tootUnknown)
			throw ParseException(petUnknownInput, std::move(retval));
	}

	void InputParser::parseFile(const std::string& fileName)
	{
		curFileName = fileName;
		valid = false;
		// Set default file names for rep and opt
		pugi::xml_document xmldoc;
		pugi::xml_parse_result result = xmldoc.load_file(fileName.c_str());
		if(!result)
		{
			std::cerr << "Error: couldn't load input file: " << fileName << ", aborting" << std::endl;
			std::cout << "Error description: " << result.description() << std::endl;
			abort();
		}
		parse(xmldoc);
		valid = true;
	}
	
	void InputParser::parseNode(const NodeInfo& ni)
	{
		curFileName = ni.getFileName();
		valid = false;
		pugi::xml_document xmldoc;
		xmldoc.load_file(curFileName.c_str());
		if(!xmldoc)
		{
			std::cerr << "Error: couldn't load input file: " << curFileName << ", aborting" << std::endl;
			abort();
		}
		// Find the specified node
		pugi::xml_node rootNode, searchNode;
		for(rootNode = xmldoc.child(ni.getTypeName().c_str()); rootNode && !searchNode; 
				rootNode = rootNode.next_sibling(ni.getTypeName().c_str()))
		{
			std::string tag = readStringAttribute(rootNode, "tag");
			if(tag == ni.getTag())
				searchNode = rootNode;
		}
		// Parse
		if(searchNode)
			parse(searchNode);
		else
		{
			std::cerr << "Error: Couldn't find node " << ni.getTypeName() << " with tag " << ni.getTag() << " in file " << curFileName << std::endl;
			abort();
		}
		valid = true;
	}
}	
}

