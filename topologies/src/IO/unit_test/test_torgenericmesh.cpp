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

#define CATCH_CONFIG_MAIN

#include "inputloaderrep.h"
#include "catch.hpp"
#include <string>

using namespace Topologies;

TEST_CASE("Testing input parsing using TORGenericMesh","[TORGenericMesh]")
{
	using namespace InputLoader;
	RepNodeInfo testRNI;
	REQUIRE(testRNI.getNodeName() == "representation");
	pugi::xml_document xmldoc;
	std::string fileName("testmesh2d.xml");
	xmldoc.load_file(fileName.c_str());
	REQUIRE(xmldoc);
	pugi::xml_node rootNode = xmldoc.child("representation");
	REQUIRE(rootNode);
	TORGenericMesh defaultParser(testRNI.getTypeName());
	SECTION("Mesh2D, no default values")
	{
		testRNI.parse(rootNode, "testmesh2d.xml");
		TORGenericMesh testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getMinDensity() == Approx(1.));
		REQUIRE(testParser.getPenalPower() == Approx(2.));
		REQUIRE(testParser.getThreshold() == Approx(3.));
		REQUIRE(testParser.getFiltRad() == defaultParser.getFiltRad());
		REQUIRE(testParser.getBetaHeavi() == defaultParser.getBetaHeavi());
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
		REQUIRE(testParser.getVMTORS().torMeshType == mffPolygon);
		REQUIRE(testParser.getFileName() == "");
		// Check polygon input
		const std::vector<std::vector<Point_2_base>>& polyVec = testParser.getPolyVec();
		REQUIRE(polyVec.size() == 2);
		REQUIRE(polyVec[0].size() == 3);
		REQUIRE(polyVec[1].size() == 4);
		REQUIRE(polyVec[0][0].x() == Approx(0.));
		REQUIRE(polyVec[0][0].y() == Approx(0.));
		REQUIRE(polyVec[0][1].x() == Approx(1.));
		REQUIRE(polyVec[0][1].y() == Approx(0.));
		REQUIRE(polyVec[0][2].x() == Approx(0.));
		REQUIRE(polyVec[0][2].y() == Approx(1.));
		REQUIRE(polyVec[1][0].x() == Approx(0.1));
		REQUIRE(polyVec[1][0].y() == Approx(0.1));
		REQUIRE(polyVec[1][1].x() == Approx(0.2));
		REQUIRE(polyVec[1][1].y() == Approx(0.1));
		REQUIRE(polyVec[1][2].x() == Approx(0.2));
		REQUIRE(polyVec[1][2].y() == Approx(0.2));
		REQUIRE(polyVec[1][3].x() == Approx(0.1));
		REQUIRE(polyVec[1][3].y() == Approx(0.2));
		// Check mesh params
		const GeometryTranslation::MesherData& meshParams = testParser.getMeshParams();
		REQUIRE(meshParams.triMeshEdgeAngle == Approx(10.));
		REQUIRE(meshParams.triMeshEdgeSize == Approx(0.1));
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("HeavisideMesh2D, with default values")
	{
		testRNI.parse(rootNode, "testmesh2d.xml");
		TORGenericMesh testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getMinDensity() == defaultParser.getMinDensity());
		REQUIRE(testParser.getThreshold() == defaultParser.getThreshold());
		REQUIRE(testParser.getPenalPower() == defaultParser.getPenalPower());
		REQUIRE(testParser.getBetaHeavi() == Approx(1.));
		REQUIRE(testParser.getFiltRad() == Approx(2.));
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
		REQUIRE(testParser.getVMTORS().torMeshType == mffExodus);
		REQUIRE(testParser.getFileName() == "test.txt");
		std::vector<std::pair<unsigned, double>> blockVec = testParser.getFixedBlockVec();
		REQUIRE(blockVec.size() == 2);
		REQUIRE(blockVec[0].first == 2);
		REQUIRE(blockVec[0].second == Approx(1.));
		REQUIRE(blockVec[1].first == 3);
		REQUIRE(blockVec[1].second == Approx(0.5));
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Mesh3d, no defaults")
	{
		testRNI.parse(rootNode, "testmesh2d.xml");
		TORGenericMesh testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getMinDensity() == defaultParser.getMinDensity());
		REQUIRE(testParser.getThreshold() == defaultParser.getThreshold());
		REQUIRE(testParser.getPenalPower() == defaultParser.getPenalPower());
		REQUIRE(testParser.getFiltRad() == defaultParser.getFiltRad());
		REQUIRE(testParser.getBetaHeavi() == defaultParser.getBetaHeavi());
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
		REQUIRE(testParser.getVMTORS().torMeshType == mffSTL);
		REQUIRE(testParser.getFileName() == "test.txt");
		// Check mesh params
		const GeometryTranslation::MesherData& meshParams = testParser.getMeshParams();
		REQUIRE(meshParams.tetMeshEdgeSize == Approx(1.));
		REQUIRE(meshParams.tetMeshEdgeSize == Approx(1.));
		REQUIRE(meshParams.tetMeshFacetAngle == Approx(2.));
		REQUIRE(meshParams.tetMeshFacetSize == Approx(3.));
		REQUIRE(meshParams.tetMeshFacetDistance == Approx(4.));
		REQUIRE(meshParams.tetMeshCellRadiusEdgeRatio == Approx(5.));
		REQUIRE(meshParams.tetMeshCellSize == Approx(6.));
	}
	rootNode = rootNode.next_sibling("representation");
	REQUIRE(rootNode);
	SECTION("Heaviside3D, with default values")
	{
		testRNI.parse(rootNode, "testmesh2d.xml");
		TORGenericMesh testParser(testRNI.getTypeName());
		testParser.parseNode(testRNI);
		REQUIRE(testParser.getMinDensity() == defaultParser.getMinDensity());
		REQUIRE(testParser.getThreshold() == defaultParser.getThreshold());
		REQUIRE(testParser.getPenalPower() == defaultParser.getPenalPower());
		REQUIRE(testParser.getFiltRad() == defaultParser.getFiltRad());
		REQUIRE(testParser.getBetaHeavi() == defaultParser.getBetaHeavi());
		REQUIRE(testParser.getRepName() == testRNI.getTypeName());
		REQUIRE(testParser.getVMTORS().torMeshType == mffExodus);
		REQUIRE(testParser.getFileName() == "test.txt");
	}
}
