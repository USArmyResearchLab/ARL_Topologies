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

#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include "tofemobjfun.h"
#include "REP/topoptrep.h"
#include "femproblem.h"
#include "REP/tomesh.h"
#include "IO/exotxtmeshloader.h"
#include "UTIL/helper.h"
#include "point2d.h"
#include "point3d.h"

using namespace Topologies;

TOFEMObjFun::TOFEMObjFun(const std::string& inpFileName)
{
	FEMInputLoader inptParser;
	inptParser.parseFile(inpFileName);
	// Convert to Lame parameters
	double E = inptParser.getYoungsMod(), nu = inptParser.getPoissonsRatio();
	double lambda = E*nu/((1. + nu)*(1. - 2.*nu)), mu = 0.5*E/(1. + nu);
	std::vector<double> matParams(3);
	matParams[0] = inptParser.getDensity();
	matParams[1] = lambda;
	matParams[2] = mu;
	baseMat = GenericMaterial(matParams);
	// Copy remaining input parameters
	maxDisplacement = inptParser.getMaxDisplacement();
	volfrac = inptParser.getVolFrac();
	dim = inptParser.getDim();
	// Load boundary conditions
	bcVec = inptParser.getBCVec();
	lcVV = inptParser.getLCVec();
}

TOFEMObjFun::~TOFEMObjFun()
{
}

std::unique_ptr<TOMesh> TOFEMObjFun::getTOMesh(const TopOptRep& inTOR) const
{
	if(dim == 2)
		return inTOR.get2DMesh();
	return inTOR.get3DVolumeMesh();
}

std::unique_ptr<FEMProblem> TOFEMObjFun::setupAndSolveFEM(const TopOptRep& inTOR, std::size_t kLoad) const
{
	// Set up mesh and fem problem
	std::unique_ptr<TOMesh> resMesh = getTOMesh(inTOR);
	assert(resMesh);
	std::unique_ptr<FEMProblem> upProb(new FEMProblem(resMesh.get(), baseMat));
	// Convert boundary conditions
	std::vector<ExoBC> exoBCVec = generateExoBCVec(resMesh.get(), kLoad);
	// Solve problem and get compliance
	upProb->changeBoundaryConditionsTo(exoBCVec);
	return upProb;
}

std::pair<double, bool> TOFEMObjFun::operator()(const TopOptRep& inTOR) const
{
	std::pair<double, bool> rval(0., true);
	for(std::size_t k = 0; k < lcVV.size() && rval.second; ++k)
	{
		// Loop over each load case and sum compliance
		std::unique_ptr<FEMProblem> upFEM = setupAndSolveFEM(inTOR, k);
		std::pair<double, bool> tmp = upFEM->computeCompliance(); 
		rval.first += tmp.first;
		rval.second &= tmp.second;
	}
	return rval;
}

void TOFEMObjFun::f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{	
	std::pair<double, bool> res = (*this)(inTOR);
	std::vector<double> ofvVec(2);
	ofvVec[0] = res.first;
	ofvVec[1] = inTOR.computeVolumeFraction();
	bool valid = res.second;
	if(res.first > maxDisplacement)
		valid = false;
	outRes = std::make_pair(ofvVec, valid);
}

void TOFEMObjFun::g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	std::vector<std::map<std::size_t, double>> dTOR = inTOR.diffRep();
	if(dTOR.empty()) // No analytical derivative is available, use finite differences
		TopOptObjFun::g(inTOR, outRes);
	else
		g(inTOR, dTOR, outRes);
}

void TOFEMObjFun::g(const TopOptRep& inTOR, const std::vector<std::map<std::size_t, double>>& dTOR,
    std::pair<std::vector<double>, bool>& outRes) const
{
	// Gradient consists of partial derivatives of inTOR multiplied by those of the FEM problem
	outRes.first = std::vector<double>(dTOR.size(), 0.);
	outRes.second = true;
	for(std::size_t k = 0; k < lcVV.size() && outRes.second; ++k)
	{
		std::unique_ptr<FEMProblem> upFEM = setupAndSolveFEM(inTOR, k);
		outRes.second &= upFEM->validRun();
		std::unique_ptr<TOMesh> resMesh = getTOMesh(inTOR);
		assert(resMesh);
		outRes.first = HelperNS::vecSum(outRes.first, upFEM->gradCompliance(resMesh.get(), dTOR));
	}
}

// Constraints
void TOFEMObjFun::c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	double val = inTOR.computeVolumeFraction() - volfrac;
	std::vector<double> ofvVec(1, val);
	outRes = std::make_pair(ofvVec, true);
}

void TOFEMObjFun::gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes) const
{
	std::vector<double> dVF = inTOR.computeGradVolumeFraction();
	if(dVF.empty())
		TopOptObjFun::gc(inTOR, outRes);
	else
		outRes = std::make_pair(std::move(dVF), true);
}

#ifdef USE_MPI
std::pair<double, bool> TOFEMObjFun::operator()(const TopOptRep& inTOR, MPI::Comm& communicator) const
{
	// Not implemented
	return (*this)(inTOR);
}

void TOFEMObjFun::f(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	f(inTOR, outRes);
}

void TOFEMObjFun::c(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	c(inTOR, outRes);
}

void TOFEMObjFun::g(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	std::vector<std::map<std::size_t, double>> dTOR = inTOR.diffRep();
	if(dTOR.empty()) // No analytical derivative is available, use finite differences
		TopOptObjFun::g(inTOR, outRes, communicator);
	else
		g(inTOR, dTOR, outRes);
}

void TOFEMObjFun::gc(const TopOptRep& inTOR, std::pair<std::vector<double>, bool>& outRes, MPI::Comm& communicator) const
{
	std::vector<double> dVF = inTOR.computeGradVolumeFraction();
	if(dVF.empty())
		TopOptObjFun::gc(inTOR, outRes, communicator);
	else
		outRes = std::make_pair(std::move(dVF), true);
}

#endif

std::vector<ExoBC> TOFEMObjFun::generateExoBCVec(const TOMesh* const inMesh, std::size_t kLoad) const
{
	// First boundary conditions
	std::vector<ExoBC> outEBC;
	assert(loadCase < lcVV.size());
	const std::vector<LoadCondition<double>>& lcVec = lcVV[kLoad];
	outEBC.reserve(bcVec.size() + lcVec.size());
	for(auto const& bc : bcVec)
	{
		outEBC.push_back(ExoBC(dim));
		ExoBC& curebc = outEBC.back();
		curebc.dim = dim;
		curebc.isSupport = true;
		std::vector<bool> supps = bc.getFixedCoords();
		curebc.xsup = supps[0];
		curebc.ysup = supps[1];
		if(dim == 3)
			curebc.zsup = supps[2];
		curebc.nodeIDVec = bc.applyBC(inMesh);
	}	
	// Now load conditions
	for(auto const& lc : lcVec)
	{
		outEBC.push_back(ExoBC(dim));
		ExoBC& curebc = outEBC.back();
		curebc.dim = dim;
		curebc.isSupport = false;
		lc.applyLC(inMesh, curebc.nodeIDVec);
		std::vector<double> loadVec = lc.getLoadVec();
		assert(loadVec.size() >= 2);
		curebc.loadVec = Point3D(loadVec[0], loadVec[1], loadVec.size() > 2 ? loadVec[2] : 0.);
		curebc.ct = lc.getCoordinateSystemType();
	}
	return outEBC;
}

// Input parsing

void FEMInputLoader::parse(const pugi::xml_document& xmldoc)
{
	parse(xmldoc.child("fem"));
}	

void FEMInputLoader::parse(const pugi::xml_node& rootNode)
{
	using namespace InputLoader;
	try
	{
		// Read required inputs
		E = readDoublePCData(rootNode, "youngs_modulus");
		nu = readDoublePCData(rootNode, "poissons_ratio");
		volfrac = readDoublePCData(rootNode, "volume_fraction");
		dim = readUnsignedPCData(rootNode, "dimension");
		if(dim != 2 && dim != 3)
			throw(ParseException(petInvalidNumericalInput, "dimension must be 2 or 3"));
	}
	catch(ParseException pe)
	{
		errorMessage(pe, true);
	}
	// Optional inputs
	try{maxDisplacement = readDoublePCData(rootNode, "max_displacement");}
	catch(ParseException pe){maxDisplacement = 1e16;}
	try{rho = readDoublePCData(rootNode, "density");}
	catch(ParseException pe){rho = 1.;}
	try{meshFilename = readAndCheckFileNamePCData(rootNode, "mesh_file");}
	catch(ParseException pe){}//Only necessary for file_bc input
	try{setFileFormat(rootNode);}
	catch(ParseException pe){}// Only necessary for file_bc input
	if(!meshFilename.empty() && inputMFF == mffExodus)
	{
		// Check mesh dimensions
		unsigned probDim = ExoTxtMeshLoader::readMeshDimension(meshFilename);
		if(probDim != dim)
		{
			ParseException pe(petInvalidNumericalInput, "Wrong number of dimensions in exodus mesh file!");
			errorMessage(pe, true);
		}
	}
	// ExoBC input
	for(pugi::xml_node cn = rootNode.child("file_bc"); cn; cn = cn.next_sibling("file_bc"))
		parseExoBC(cn);
	// GeoBC input
	for(pugi::xml_node cn = rootNode.child("geo_bc"); cn; cn = cn.next_sibling("geo_bc"))
		parseGeoBC(cn);
	// Parse load cases (each load case will run a separate FEM solve)
	for(pugi::xml_node cn1 = rootNode.child("load_case"); cn1; cn1 = cn1.next_sibling("load_case"))
	{
		lcVV.push_back(std::vector<LoadCondition<double>>());
		// ExoLC input
		for(pugi::xml_node cn2 = cn1.child("file_lc"); cn2; cn2 = cn1.next_sibling("file_lc"))
			parseExoLC(cn2);
		// GeoBC input
		for(pugi::xml_node cn2 = cn1.child("geo_lc"); cn2; cn2 = cn2.next_sibling("geo_lc"))	
			parseGeoLC(cn2);	
	}
}

void FEMInputLoader::parseExoBC(const pugi::xml_node& rootNode)
{
	using namespace InputLoader;
	if(meshFilename.empty())
	{
		ParseException pe(petPCDataNotFound, "No exodus mesh input, though using exodus boundary conditions");
		errorMessage(pe, true);
	}
	bool xsup, ysup, zsup = false;
	unsigned nodeSetID;
	try
	{
		nodeSetID = readUnsignedPCData(rootNode, "node_set_id");
		xsup = readBoolPCData(rootNode, "x_support");
		ysup = readBoolPCData(rootNode, "y_support");
		if(dim == 3)
			zsup = readBoolPCData(rootNode, "z_support");
	}
	catch(ParseException pe)
	{
		errorMessage(pe, true);
	}
	bcVec.push_back(BoundaryCondition(bctUnknown, xsup, ysup, zsup, nodeSetID, inputMFF, meshFilename, dim));
}

void FEMInputLoader::parseExoLC(const pugi::xml_node& rootNode)
{
	using namespace InputLoader;
	if(meshFilename.empty())
	{
		ParseException pe(petPCDataNotFound, "No exodus mesh input, though using exodus boundary conditions");
		errorMessage(pe, true);
	}
	unsigned nodeSetID;
	std::vector<double> vecRes;
  try
  {
		nodeSetID = readUnsignedPCData(rootNode, "node_set_id");
		vecRes = readDoubleVecPCData(rootNode, "load_vector");
		if(vecRes.size() != dim)
			throw ParseException(petInvalidNumericalInput, "Wrong number of dimensions in load_vector!");
	}
  catch(ParseException pe)
  {
    errorMessage(pe, true);
  }
	CoordinateSystem::Type ct = readCoordinateSystemAttribute(rootNode.child("load_vector"));
	lcVV.back().push_back(LoadCondition<double>(bctUnknown, vecRes, nodeSetID, inputMFF, meshFilename, dim, ct));
}

CoordinateSystem::Type FEMInputLoader::parseCSType(const std::string& inCSTStr)
{
	using namespace InputLoader;
	if(inCSTStr == "cartesian")
		return CoordinateSystem::Type::cartesian;
	else if(inCSTStr == "cylindrical")
		return CoordinateSystem::Type::cylindrical;
	else if(inCSTStr == "spherical")
		return CoordinateSystem::Type::spherical;
	throw ParseException(petUnknownInput, "Unknown coordinate system type");
}

BCType FEMInputLoader::parseBCType(const std::string& inBCStr)
{
	if(inBCStr == "v_line")
		return bctVLine;
	else if(inBCStr == "h_line")
		return bctHLine;
	else if(inBCStr == "line_segment_2")
		return bctLineSeg;
	else if(inBCStr == "point_2")
		return bctPoint2D;
	else if(inBCStr == "polygon_2")
		return bctPolygon2D;
	else if(inBCStr == "xy_plane")
		return bctXYPlane;
	else if(inBCStr == "yz_plane")
		return bctYZPlane;
	else if(inBCStr == "xz_plane")
		return bctXZPlane;
	else if(inBCStr == "point_3")
		return bctPoint3D;
	else if(inBCStr == "line_segment_3")
		return bctLineSeg3D;
	return bctUnknown;
}

void FEMInputLoader::parseGeoBC(const pugi::xml_node& rootNode)
{
	using namespace InputLoader;
	BCType outBCT;
	bool xsup, ysup, zsup = false;
	std::unique_ptr<GeometricEntity> upGE;
	try
	{
		outBCT = parseBCType(rootNode.child("geometry").first_child().name());
		upGE = GeometryEntityFactory::createGeometricEntity(outBCT, rootNode.child("geometry").first_child());
		xsup = readBoolPCData(rootNode, "x_support");
		ysup = readBoolPCData(rootNode, "y_support");
		if(dim == 3)
			zsup = readBoolPCData(rootNode, "z_support");
	}
	catch(ParseException pe)
	{
		errorMessage(pe, true);
	}
	bcVec.push_back(BoundaryCondition(outBCT, xsup, ysup, zsup, std::move(upGE)));
}

CoordinateSystem::Type FEMInputLoader::readCoordinateSystemAttribute(const pugi::xml_node& rootNode) const
{
	using namespace InputLoader;
	CoordinateSystem::Type ct;
	std::string cs;
	try
	{
		cs = readStringAttribute(rootNode, "coordinate_type");
		ct = parseCSType(cs);
	}
	catch(ParseException pe)
	{
		ct = CoordinateSystem::Type::cartesian;
		if(!cs.empty())
			std::cerr << "Warning: Unknown coordinate system in load_vector: " << cs << "\nDefaulting to cartesian" << std::endl;
	}
	return ct;
}

void FEMInputLoader::parseGeoLC(const pugi::xml_node& rootNode)
{
	using namespace InputLoader;
	BCType outBCT;
	std::vector<double> vecRes;
	std::unique_ptr<GeometricEntity> upGE;
  try
  {
		outBCT = parseBCType(rootNode.child("geometry").first_child().name());
    upGE = GeometryEntityFactory::createGeometricEntity(outBCT, rootNode.child("geometry").first_child());
		vecRes = readDoubleVecPCData(rootNode, "load_vector");
    if(vecRes.size() != dim)
      throw ParseException(petInvalidNumericalInput, "Wrong number of dimensions in load_vector!");
  }
  catch(ParseException pe)
  {
    errorMessage(pe, true);
  }
	CoordinateSystem::Type ct = readCoordinateSystemAttribute(rootNode.child("load_vector"));
  lcVV.back().push_back(LoadCondition<double>(outBCT, vecRes, std::move(upGE), ct));
}

void FEMInputLoader::setFileFormat(const pugi::xml_node& rootNode)
{
	using namespace InputLoader;
	std::string retval = readStringPCData(rootNode, "file_format");
	inputMFF = parseMeshFileFormat(retval);
	if(inputMFF != mffExodus && inputMFF != mffGMSH)
		throw ParseException(petUnknownInput, std::move(retval));
}

