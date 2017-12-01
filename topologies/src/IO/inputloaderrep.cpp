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

#include "inputloaderrep.h"

namespace Topologies{
namespace InputLoader
{
	void TORGenericVolume::parse(const pugi::xml_document& xmldoc)
	{
		parse(xmldoc.child(repName.c_str()));
	}

	void TORGenericVolume::parse(const pugi::xml_node& rootNode)
	{
		// Initialize VoLMeshTORSpec
		myVMTORS = VolMeshTORSpecification(parseTORType(rootNode.name()));
		// Required inputs
		try
		{
			discSizes = readUnsignedVecPCData(rootNode, "discretization_size");
			dimensions = readDoubleVecPCData(rootNode, "dimensions");
			if(discSizes.size() != dimensions.size())
			{
				std::string errorStr("discretization_size doesn't match dimensions");
				throw ParseException(petInvalidNumericalInput, std::move(errorStr));
			}
			setMeshElementType(rootNode);
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
		// Optional inputs
		try{threshold = readDoublePCData(rootNode, "threshold");}
		catch(ParseException pe){}
		try{penalPower = readDoublePCData(rootNode, "penalty_power");}
		catch(ParseException pe){}
		try{filtRad = readDoublePCData(rootNode, "filter_radius");}
		catch(ParseException pe){}
		try{betaHeavi = readDoublePCData(rootNode, "heaviside_beta");}
		catch(ParseException pe){}
		try{minDensity = readDoublePCData(rootNode, "minimum_density");}
		catch(ParseException pe){}
		try{rank = readUnsignedPCData(rootNode, "rank");}
		catch(ParseException pe){}
	}

	void TORGenericVolume::setMeshElementType(const pugi::xml_node& rootNode)
	{
		std::string retval = readStringPCData(rootNode, "mesh_element_type");
		myMET = parseMeshElementType(retval);
		if(myMET == metUnknown)
			throw ParseException(petUnknownInput, std::move(retval));
	}

	void TORCSGTree::parse(const pugi::xml_document& xmldoc)
	{
		parse(xmldoc.child(repName.c_str()));
	}

	void TORCSGTree::parse(const pugi::xml_node& rootNode)
	{
		try
		{
			regionDimensions = readDoubleVecPCData(rootNode, "region_dimensions");
			meshSizes = readUnsignedVecPCData(rootNode, "mesh_sizes");
			shapeNums = readUnsignedVecPCData(rootNode, "num_shapes");
			nPointsPerShape = readUnsignedPCData(rootNode, "num_points_per_shape");
			if(regionDimensions.size() != meshSizes.size() || regionDimensions.size() != shapeNums.size() || 
				meshSizes.size() != shapeNums.size())
			{
				std::string errorStr("Mismatch in input parameter dimensions");
				throw ParseException(petInvalidNumericalInput, std::move(errorStr));
			}
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
		// Optional
		try{minDensity = readDoublePCData(rootNode, "minimum_density"); }
		catch(ParseException pe){}
		try{useAffine = readBoolPCData(rootNode, "use_affine");}
		catch(ParseException pe){}
		try{shapesAreHoles = readBoolPCData(rootNode, "shapes_are_holes");}
		catch(ParseException pe){}
	}

	void TORGenericMesh::parse(const pugi::xml_document& xmldoc)
	{
		parse(xmldoc.child(repName.c_str()));
	}

	void TORGenericMesh::parse(const pugi::xml_node& rootNode)
	{
		// Initialize VoLMeshTORSpec
		myVMTORS = VolMeshTORSpecification(parseTORType(rootNode.name()));
		// Required inputs
		try
		{
			setFileFormat(rootNode);
			if(myVMTORS.torMeshType != mffPolygon)
				fileName = readAndCheckFileNamePCData(rootNode, "file_name");
			else
				setPolygonList(rootNode.child("polygon_list"));
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
		// Optional
		try{minDensity = readDoublePCData(rootNode, "minimum_density"); }
		catch(ParseException pe){}
		try{threshold = readDoublePCData(rootNode, "threshold");}
		catch(ParseException pe){}
		try{penalPower = readDoublePCData(rootNode, "penalty_power");}
		catch(ParseException pe){}
		try{filtRad = readDoublePCData(rootNode, "filter_radius");}
		catch(ParseException pe){}
		try{betaHeavi = readDoublePCData(rootNode, "heaviside_beta");}
		catch(ParseException pe){}

		//Mesh parameters
		if(myVMTORS.torMeshType == mffPolygon || myVMTORS.torMeshType == mffSTL)
			meshParams.parse(rootNode.child("mesh_parameters"));
		else
		{
			// Fixed blocks
			for(pugi::xml_node cn = rootNode.child("fixed_block"); cn; cn = cn.next_sibling("fixed_block"))
			{
				std::pair<unsigned,double> fixedBlockDef;
				try{
					fixedBlockDef.first = readUnsignedAttribute(cn, "block_id");
					fixedBlockDef.second = readDoubleAttribute(cn, "value");
				}
				catch(ParseException pe) {errorMessage(pe, true);}
				fixedBlockVec.push_back(fixedBlockDef);
			}
		}
	}

	void TORGenericMesh::setFileFormat(const pugi::xml_node& rootNode)
	{
		std::string retval = readStringPCData(rootNode, "file_format");
		myVMTORS.torMeshType = parseMeshFileFormat(retval);
		if(myVMTORS.torMeshType == mffUnknown)
			throw ParseException(petUnknownInput, std::move(retval));
	}

	void TORGenericMesh::setPolygonList(const pugi::xml_node& rootNode)
	{
		for(pugi::xml_node cn = rootNode.child("polygon_2"); cn; cn = cn.next_sibling("polygon_2"))
			ptVecs.push_back(parsePointVec(cn));
	}
		
	std::vector<Point_2_base> TORGenericMesh::parsePointVec(const pugi::xml_node& rootNode)
	{
		std::vector<Point_2_base> ptVec;
		for(pugi::xml_node cn = rootNode.child("point_2"); cn; cn = cn.next_sibling("point_2"))
		{
			double x = readDoubleAttribute(cn, "x");
			double y = readDoubleAttribute(cn, "y");
			ptVec.push_back(Point_2_base(x, y));
		}
		return ptVec;
	}

	void CGALMesh::parse(const pugi::xml_document& xmldoc)
	{
		parse(xmldoc.child("mesh_parameters"));
	}

	void CGALMesh::parse(const pugi::xml_node& rootNode)
	{
		try
		{
			dim = readUnsignedAttribute(rootNode, "dim");
			if(dim == 2)
			{
				inputParams.triMeshEdgeSize = readDoublePCData(rootNode, "tri_mesh_edge_size");
				inputParams.triMeshEdgeAngle = readDoublePCData(rootNode, "tri_mesh_edge_angle");
			}
			else if(dim == 3)
			{
				inputParams.tetMeshEdgeSize = readDoublePCData(rootNode, "tet_mesh_edge_size");
				inputParams.tetMeshFacetAngle = readDoublePCData(rootNode, "tet_mesh_facet_angle");
				inputParams.tetMeshFacetSize = readDoublePCData(rootNode, "tet_mesh_facet_size");
				inputParams.tetMeshFacetDistance = readDoublePCData(rootNode, "tet_mesh_facet_distance");
				inputParams.tetMeshCellRadiusEdgeRatio = readDoublePCData(rootNode, "tet_mesh_cell_radius_edge_ratio");
				inputParams.tetMeshCellSize = readDoublePCData(rootNode, "tet_mesh_cell_size");
			}
			else
			{
				std::string errorMsg("Wrong dimension read in mesh_parameters, must be 2 or 3");
				throw ParseException(petInvalidNumericalInput, std::move(errorMsg));
			}
		}
		catch(ParseException pe)
		{
			errorMessage(pe, true);
		}
	}
}
}

