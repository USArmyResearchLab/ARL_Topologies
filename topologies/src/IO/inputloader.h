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

#ifndef INPUTLOADER_H
#define INPUTLOADER_H

#include <string>
#include <vector>
#include <pugixml.hpp>
#include "geometrytranslation.h"

namespace Topologies{
class TopOptRep;

//! Enumeration for different types of initial guesses
enum InitialGuessType{igtConstant, igtConstantWithNoise, igtRandom, igtFile, igtUnknown};
//! Enumeration for output types
enum OutputType{otSurface, otVolume, otExtrude, otObjFunRes, otRawData, otUnknown};
//! Enumeration for output file formats
enum OutputFileFormat{offMatlab, offSTL, offVTK, offGMSH, offDefault};
//! Enumeration for mesh element types
enum MeshElementType {metTri, metQuad, metTet, metHex, metUnknown};
//! Enumeration for mesh file formats
enum MeshFileFormat{mffExodus, mffGMSH, mffPolygon, mffSTL, mffStructured, mffUnknown};
//! Enumeration to define 2d or 3d
enum DimensionType {dt2d, dt3d};
//! Enumeration for penalization methods
enum PenalizationType {ptSIMP, ptRAMP};
//! Enumeration for projection methods
enum ProjectionType {ptNone, ptThresholdHeavi};
//! Enumeration for filter types
enum FilterType {ftNone, ftLinear};
//! Enumeration for unknown location (nodal vs. element)
enum UnknownLocation {ulElement, ulNode};
//! Struct to define options for volumetric mesh-based topology representations
struct VolMeshTORSpecification
{
	//! Default constructor, defines tortPixel
	VolMeshTORSpecification(); 
	//! Constructor taking TORType, defines what information it can
	VolMeshTORSpecification(TORType inTORT);
	//! Constructor setting all values
	VolMeshTORSpecification(DimensionType inDT, MeshFileFormat inMFF, UnknownLocation inUL, FilterType inFT, 
		PenalizationType inPenType, ProjectionType inProjType);
	//! Constructor using a vector of ints, size must match number of enums
	explicit VolMeshTORSpecification(const std::vector<int>& params);
	//! Returns a vector of ints, cast from the enums
	std::vector<int> toVec() const;
	//! Returns the number of options, using the result of the function toVec
	std::size_t size() const;
	// Options
	DimensionType torDimType;
	MeshFileFormat torMeshType;
	UnknownLocation torUnknownLocation;
	FilterType torFilterType;
	PenalizationType torPenalizationType;
	ProjectionType torProjectionType;
};

//! A collection of functions and namespaces for loading input
/*! This namespace handles all file input for Topologies.  The name space contains several
*  nested namespaces for loading different aspects of the method, such as TopOptReps,
*  optimizers, etc.  Each sub-namespace contains variables to hold inputs and 
*  functions that parse the input and interface with the Parser class.
*/
namespace InputLoader
{
	//! Loads a TopOptRep from a data file
	/*! This function will read a TOR data file that was output by Topologies.
	*  This can be used to load a save file if an optimization run was terminated
	*  prematurely or it can be used to perform other output functions on a previously
	*  finished run.
	*/
	void inputTORData(TopOptRep* const outTOR, const std::string& fileName);
	//! Convert string to InitialGuessType
	InitialGuessType parserIGType(const std::string& strTOOType);
	//! Convert string to TORType
	TORType parseTORType(const std::string& strTORType);
	//! Convert string to TOOType
	TOOType parseTOOType(const std::string& strTOOType);
	//! Convert string to OutputType
	OutputType parseOutputType(const std::string& strOType);
	//! Convert string to OutputFileFormat
	OutputFileFormat parseOutputFileFormat(const std::string& strOType);
	//! Convert string to MeshElementType
	MeshElementType parseMeshElementType(const std::string& strOType);
	//! Convert string to MeshElementType
	MeshFileFormat parseMeshFileFormat(const std::string& strOType);
	//! Parsing exception types
	enum ParseExceptionType{petPCDataNotFound, petAttributeNotFound, petWrongNodeType, petUnknownInput, petInvalidNumericalInput};
	//! Struct to output meaningful information in the event of an input parsing error
	struct ParseException
	{
		ParseException(const ParseExceptionType& inPET, std::string&& inInfo) :
			pet(inPET), info(std::move(inInfo)) {}
		ParseExceptionType pet;
		std::string info;
	};

	//! Returns name of ParseException for output purposes
	std::string getExceptionName(const ParseException& inPE);
	//! Output error message
	void errorMessage(const ParseException& inPE, bool willExit);

	//! Read string from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	std::string readStringPCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path);
	//! Read string from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when attribute is not found
	std::string readStringAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName);
	//! Read double from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	double readDoublePCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path);
	//! Read double from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	double readDoubleAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName);
	//! Read unsigned int from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	unsigned readUnsignedPCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path);
	//! Read unsigned from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	unsigned readUnsignedAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName);
	//! Read int from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	int readIntPCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path);
	//! Read int from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	int readIntAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName);
	//! Read bool from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	bool readBoolPCData(const pugi::xml_node& rootNode, const std::vector<std::string>& path);
	//! Read bool from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	bool readBoolAttribute(const pugi::xml_node& rootNode, const std::vector<std::string>& path, const std::string& attrName);
	//! Read string from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	std::string readStringPCData(const pugi::xml_node& rootNode, const std::string& path);
	//! Read string from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	std::string readStringAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName);
	//! Read double from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	double readDoublePCData(const pugi::xml_node& rootNode, const std::string& path);
	//! Read double from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	double readDoubleAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName);
	//! Read unsigned int from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	unsigned readUnsignedPCData(const pugi::xml_node& rootNode, const std::string& path);
	//! Read unsigned from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	unsigned readUnsignedAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName);
	//! Read int from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	int readIntPCData(const pugi::xml_node& rootNode, const std::string& path);
	//! Read int from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	int readIntAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName);
	//! Read int from specified path in XML tree
	//! @throws ParseException Throws exception when path is not found
	bool readBoolPCData(const pugi::xml_node& rootNode, const std::string& path);
	//! Read int from specified attribute at path in XML tree
	//! @throws ParseException Throws exception when path is not found
	bool readBoolAttribute(const pugi::xml_node& rootNode, const std::string& path, const std::string& attrName);
	//! Read string from specified attribute in XML tree
	//! @throws ParseException Throws exception when path is not found
	std::string readStringAttribute(const pugi::xml_node& rootNode, const std::string& attrName);
	//! Read double from specified attribute in XML tree
	//! @throws ParseException Throws exception when path is not found
	double readDoubleAttribute(const pugi::xml_node& rootNode, const std::string& attrName);
	//! Read unsigned from specified attribute in XML tree
	//! @throws ParseException Throws exception when path is not found
	unsigned readUnsignedAttribute(const pugi::xml_node& rootNode, const std::string& attrName);
	//! Read int from specified attribute in XML tree
	//! @throws ParseException Throws exception when path is not found
	int readIntAttribute(const pugi::xml_node& rootNode, const std::string& attrName);
	//! Read bool from specified attribute in XML tree
	//! @throws ParseException Throws exception when path is not found
	bool readBoolAttribute(const pugi::xml_node& rootNode, const std::string& attrName);
	//! Read string from XML tree
	//! @throws ParseException Throws exception when path is not found
	std::string readStringPCData(const pugi::xml_node& rootNode);
	//! Read double from XML tree
	//! @throws ParseException Throws exception when path is not found
	double readDoublePCData(const pugi::xml_node& rootNode);
	//! Read unsigned int XML tree
	//! @throws ParseException Throws exception when path is not found
	unsigned readUnsignedPCData(const pugi::xml_node& rootNode);
	//! Read int from XML tree
	//! @throws ParseException Throws exception when path is not found
	int readIntPCData(const pugi::xml_node& rootNode);
	//! Read bool from XML tree
	//! @throws ParseException Throws exception when path is not found
	bool readBoolPCData(const pugi::xml_node& rootNode);
	//! Read vector of double from the XML tree
	//! @throws ParseException Throws exception when path is not found
	std::vector<double> readDoubleVecPCData(const pugi::xml_node& rootNode, const std::string& childName);
	//! Read vector of unsigned int from XML tree
	//! @throws ParseException Throws exception when path is not found
	std::vector<unsigned> readUnsignedVecPCData(const pugi::xml_node& rootNode, const std::string& childName);
	//! Read vector of int from XML tree
	//! @throws ParseException Throws exception when path is not found
	std::vector<int> readIntVecPCData(const pugi::xml_node& rootNode, const std::string& childName);
	//! Returns the node path to the current node
	std::vector<std::string> getNodePath(const pugi::xml_node& rootNode);

	//! Node info class holds information for nodes that need to be searched later
	class NodeInfo
	{
	public:
		virtual ~NodeInfo() {}
		//! Returns the name of the node in the derived class
		virtual std::string getNodeName() const = 0;
		//! Parse an XML node
		//! @throws ParseExcxeption Throws when node is not valid and/or attributes not found
		void parse(const pugi::xml_node& rootNode, const std::string& curFileName);
		//! Parses the node type in derived classes
		//! @throws ParseExcxeption Throws when node is not valid and/or attributes not found
		virtual void parseType(const pugi::xml_node& rootNode) = 0;

		//! Returns the path to the node
		const std::vector<std::string>& getPath() const {return path;}
		//! Returns the file name containing this node
		const std::string& getFileName() const {return fileName;}
		//! Returns the node's tag
		const std::string& getTag() const {return tag;}
		//! Returns the node's name
		const std::string& getTypeName() const {return name;}
	protected:
		std::vector<std::string> path;
		std::string tag, fileName, name;
	};

	//! Derived class OptNodeInfo holds information for a TopOptRep node
	class OptNodeInfo : public NodeInfo
	{
	public:
		OptNodeInfo() : nodeTOOT(tootUnknown) {}
		virtual ~OptNodeInfo() {}
		virtual void parseType(const pugi::xml_node& rootNode);
		virtual std::string getNodeName() const {return "optimizer";}

		//! Returns the TOOType for this node
		TOOType getType() const {return nodeTOOT;}
	private:
		TOOType nodeTOOT;
	};

	//! Derived class RepNodeInfo holds information for a TopOptRep node
	class RepNodeInfo : public NodeInfo
	{
	public:
		RepNodeInfo(TORType inTORT, const NodeInfo& niData) : NodeInfo(niData), nodeTORT(inTORT) {}
		RepNodeInfo() : nodeTORT(tortUnknown) {}
		virtual ~RepNodeInfo() {}
		virtual void parseType(const pugi::xml_node& rootNode);
		virtual std::string getNodeName() const {return "representation";}

		//! Returns the TORType for this node
		TORType getType() const {return nodeTORT;}
	private:
		TORType nodeTORT;
	};

	//! Abstract base class that acts as an interface with pugixml
	class InputParser
	{
	public:
		InputParser() : valid(false) {}
		virtual ~InputParser(){}
		//! Parse the xml document, loading all input options
		virtual void parse(const pugi::xml_document& xmldoc) = 0;
		//! Parse the xml document starting at rootNode, loading all input options
		virtual void parse(const pugi::xml_node& rootNode) = 0;
		//! Parse the xml document, loading all input options
		void parseFile(const std::string& fileName);
		//! Parse the xml document, starting from the specified node
		void parseNode(const NodeInfo& ni);
		//! Returns whether or not an input file has been parsed successfully
		virtual bool isValid() const {return *this;}
		//! Returns whether or not an input file has been parsed successfully
		virtual operator bool() const {return valid;}
	protected:
		bool valid;
		std::string curFileName;
	};
}
}
#endif

