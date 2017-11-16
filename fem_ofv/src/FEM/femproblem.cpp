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

#include <Eigen/IterativeLinearSolvers>
//#include <Eigen/SparseLU>
#include "femproblem.h"
#include "element.h"
#include "lintri.h"
#include "linquad.h"
#include "lintetra.h"
#include "trilinhex.h"
#include "mesh2d.h"
#include "mesh3d.h"
#include "REP/tomesh.h"

FEMProblem::FEMProblem(const Topologies::TOMesh* const inMesh, const Topologies::GenericMaterial& baseMat) :
	numFreeDOFs(0),
	invalid(false)
{
	assert(inMesh);
	assert(baseMat.getNumParameters() > 2); // Assumes 3 material properties: density and 2 Lame parameters
	dim = inMesh->dimNum();
	if(dim == 2)
		probMesh = std::unique_ptr<FEMMesh>(new Mesh2D(inMesh, baseMat));
	else if(dim == 3)
		probMesh = std::unique_ptr<FEMMesh>(new Mesh3D(inMesh, baseMat));
}

FEMProblem::~FEMProblem()
{
}

FEMProblem::FEMProblem(const FEMProblem& copy) :
	probMesh(std::unique_ptr<FEMMesh>(copy.probMesh->clone())),
	pFEMMatrix(copy.pFEMMatrix ? std::unique_ptr<EigenSparseMat>(new EigenSparseMat(*copy.pFEMMatrix)) : nullptr),
	fixedDOFs(copy.fixedDOFs),
	bfRemapVec(copy.bfRemapVec),
	loadVec(copy.loadVec),
	pVVec(copy.pVVec ? std::unique_ptr<EigenVector>(new EigenVector(*copy.pVVec)) : nullptr),
	pForce(copy.pForce ? std::unique_ptr<EigenVector>(new EigenVector(*copy.pForce)) : nullptr),
	dim(copy.dim),
	numFreeDOFs(copy.numFreeDOFs),
	invalid(copy.invalid)
{
}

FEMProblem::FEMProblem(FEMProblem&& copy)
{
	swap(copy);
}

FEMProblem FEMProblem::operator=(FEMProblem copy)
{
	swap(copy);
	return *this;
}

void FEMProblem::swap(FEMProblem& arg2)
{
	probMesh.swap(arg2.probMesh);
	pFEMMatrix.swap(arg2.pFEMMatrix);
	pVVec.swap(arg2.pVVec);
	pForce.swap(arg2.pForce);
	std::swap(fixedDOFs, arg2.fixedDOFs);
	std::swap(bfRemapVec, arg2.bfRemapVec);
	std::swap(loadVec, arg2.loadVec);
	std::swap(dim, arg2.dim);
	std::swap(numFreeDOFs, arg2.numFreeDOFs);
	std::swap(invalid, arg2.invalid);
}

Point3D FEMProblem::getMeshPoint(std::size_t kn) const
{
	if(probMesh->getDim() == 2)
	{
		Point2D pt = *dynamic_cast<Mesh2D*>(probMesh.get())->getNode(kn);
		return Point3D(pt.x, pt.y, 0.);
	}
	return *dynamic_cast<Mesh3D*>(probMesh.get())->getNode(kn);
}

void FEMProblem::changeBoundaryConditionsTo(const std::vector<ExoBC>& bcVec)
{
	std::size_t nnodes = probMesh->getNumUnknowns();
	fixedDOFs = std::vector<bool>(dim*nnodes, false);
	loadVec.clear();
	for(auto const& curBC : bcVec)
	{
		for(auto kn : curBC.nodeIDVec)
		{
			// Set up supports
			if(curBC.isSupport)
			{
				fixedDOFs[kn] = curBC.xsup;
				fixedDOFs[kn + nnodes] = curBC.ysup;
				if(dim == 3)
					fixedDOFs[kn + 2*nnodes] = curBC.zsup;
			}
			else // Set up loads
			{
				Point3D cartVec = CoordinateSystem::convertVector(curBC.loadVec, getMeshPoint(kn), curBC.ct);
				loadVec[kn] = cartVec.x;
				loadVec[kn + nnodes] = cartVec.y;
				if(dim == 3)
					loadVec[kn + 2*nnodes] = cartVec.z;
			}
		}
	}
	// Save global row numbers for compressed system (without fixed dofs)
	bfRemapVec.resize(dim*nnodes);
	std::size_t curBF = 0;
	for(std::size_t k = 0; k < fixedDOFs.size(); ++k)
	{
		if(!fixedDOFs[k])
			bfRemapVec[k] = curBF++;
	}
	numFreeDOFs = curBF;
	solveProblem();
}

void FEMProblem::solveProblem()
{
	// Set up vectors and matrices
	setMatrix();
	setVector();
	// Matrix stuff:
	// Save RHS to use for computation of compliance
	pForce = std::unique_ptr<EigenVector>(new EigenVector(*pVVec));
	invalid = false;
	// Solve
//	Eigen::BiCGSTAB<EigenSparseMat> solver;
	Eigen::ConjugateGradient<EigenSparseMat> solver;
	solver.compute(*pFEMMatrix);
	solver.setTolerance(1e-14);
	unsigned niters = 10000;
	solver.setMaxIterations(niters);
	*pVVec = solver.solve(*pForce);
/* Direct solver
	Eigen::SparseLU<EigenSparseMat> solver;
	pFEMMatrix->makeCompressed();
	solver.analyzePattern(*pFEMMatrix);
	solver.factorize(*pFEMMatrix);
	*pVVec = solver.solve(*pForce);*/

	if(solver.iterations() >= niters)
		invalid = true;
}

void FEMProblem::setMatrix()
{
	std::size_t numNodes = probMesh->getNumUnknowns();
	std::vector<EigenT> sparseVec;
	sparseVec.reserve(10*numFreeDOFs); //Size of this depends on matrix bandwidth, which is unknown as it depends on mesh connectivity
  std::size_t nelems = probMesh->getNumElements();
  for(std::size_t ielem = 0; ielem < nelems; ++ielem)
  {
    Eigen::MatrixXd elemMat = probMesh->getElementMatrix(ielem);
    assembleMatrix(sparseVec, ielem, elemMat, numNodes);
  }
	pFEMMatrix = std::unique_ptr<EigenSparseMat>(new EigenSparseMat(numFreeDOFs,numFreeDOFs));
	pFEMMatrix->setFromTriplets(sparseVec.begin(), sparseVec.end());
}

void FEMProblem::assembleMatrix(std::vector<EigenT>& rseMat, const std::size_t kelem, const Eigen::MatrixXd& elemMat, std::size_t numUnk) const
{
	if(dim == 2)
		assembleMatrix2D(rseMat, kelem, elemMat, numUnk);
	else if(dim == 3)
		assembleMatrix3D(rseMat, kelem, elemMat, numUnk);
}

void FEMProblem::assembleMatrix2D(std::vector<EigenT>& rseMat, const std::size_t kelem, const Eigen::MatrixXd& elemMat, std::size_t numUnk) const
{
	const unsigned numElemNodes = probMesh->getNumElementNodes(kelem);
	for(unsigned iTst = 0; iTst < numElemNodes; iTst++)
	{
		std::size_t gTst = probMesh->getGlobalBF(kelem, iTst);
		for(unsigned iBas = 0; iBas < numElemNodes; iBas++)
		{
			std::size_t gBas = probMesh->getGlobalBF(kelem, iBas);
			for(unsigned tstdof = 0; tstdof < dim; ++tstdof)
			{
				for(unsigned basdof = 0; basdof < dim; ++basdof)
				{
					bool isFree = !(fixedDOFs[gTst + tstdof*numUnk] || fixedDOFs[gBas + basdof*numUnk]);
					if(fabs(elemMat(iTst + tstdof*numElemNodes, iBas + basdof*numElemNodes)) > 1.e-16 && isFree)
					{
						rseMat.push_back(EigenT(bfRemapVec[gTst + tstdof*numUnk],
																		bfRemapVec[gBas + basdof*numUnk],
																		elemMat(iTst + tstdof*numElemNodes, iBas + basdof*numElemNodes)));
					}
				}
			}
		}
	}
}

void FEMProblem::assembleMatrix3D(std::vector<EigenT>& rseMat, const std::size_t kelem, const Eigen::MatrixXd& elemMat, std::size_t numUnk) const
{
	const unsigned numElemNodes = probMesh->getNumElementNodes(kelem);
	for(unsigned iTst = 0; iTst < numElemNodes; iTst++)
	{
		std::size_t gTst = probMesh->getGlobalBF(kelem, iTst);
		for(unsigned iBas = 0; iBas < numElemNodes; iBas++)
		{
			std::size_t gBas = probMesh->getGlobalBF(kelem, iBas);
			for(unsigned tstOff = 0; tstOff < dim; ++tstOff)
			{
				for(unsigned basOff = 0; basOff < dim; ++basOff)
				{
					bool isFree = !(fixedDOFs[gTst + tstOff*numUnk] || fixedDOFs[gBas + basOff*numUnk]);
					if(fabs(elemMat(dim*iTst + tstOff, dim*iBas + basOff)) > 1.e-16 && isFree)
					{
						rseMat.push_back(EigenT(bfRemapVec[gTst + tstOff*numUnk],
																		bfRemapVec[gBas + basOff*numUnk],
																		elemMat(dim*iTst + tstOff, dim*iBas + basOff)));
					}
				}
			}
		}
	}
}

void FEMProblem::setVector()
{
  pVVec = std::unique_ptr<EigenVector>(new EigenVector(numFreeDOFs));
	pVVec->setZero();
	EigenVector& rVVec = *pVVec;
	for(std::map<std::size_t,double>::const_iterator cit = loadVec.begin(); cit != loadVec.end(); ++cit)
		rVVec(bfRemapVec[cit->first]) = cit->second;
}

std::pair<double, bool> FEMProblem::computeCompliance()
{
	if(invalid)
	{
		std::cout << "Warning, FEM OFV not valid" << std::endl;
		return std::pair<double, bool>(1e6, false);
	}
	double sum = pForce->adjoint()*(*pVVec);
	return std::pair<double, bool>(sum, true);
}

double FEMProblem::elementCompliance(const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const
{
	if(dim == 2)
		return elementCompliance2D(kelem, elemMat, numUnk);
	return elementCompliance3D(kelem, elemMat, numUnk);
}

double FEMProblem::elementCompliance2D(const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const
{
	double res = 0.;
	const unsigned numElemNodes = probMesh->getNumElementNodes(kelem);
	for(unsigned iTst = 0; iTst < numElemNodes; iTst++)
	{
		std::size_t gTst = probMesh->getGlobalBF(kelem, iTst);
		for(unsigned iBas = 0; iBas < numElemNodes; iBas++)
		{
			std::size_t gBas = probMesh->getGlobalBF(kelem, iBas);
			for(unsigned tstdof = 0; tstdof < dim; ++tstdof)
			{
				std::size_t m = gTst + tstdof*numUnk;
				for(unsigned basdof = 0; basdof < dim; ++basdof)
				{
					std::size_t n = gBas + basdof*numUnk;
					if(!(fixedDOFs[m] || fixedDOFs[n]))
						res += elemMat(iTst + tstdof*numElemNodes, iBas + basdof*numElemNodes)*(*pVVec)(bfRemapVec[m])*(*pVVec)(bfRemapVec[n]);
				}
			}
		}
	}
	return res;
}

double FEMProblem::elementCompliance3D(const std::size_t kelem, const EigenDenseMat& elemMat, std::size_t numUnk) const
{
	double res = 0.;
	const unsigned numElemNodes = probMesh->getNumElementNodes(kelem);
	for(unsigned iTst = 0; iTst < numElemNodes; iTst++)
	{
		std::size_t gTst = probMesh->getGlobalBF(kelem, iTst);
		for(unsigned iBas = 0; iBas < numElemNodes; iBas++)
		{
			std::size_t gBas = probMesh->getGlobalBF(kelem, iBas);
			for(unsigned tstOff = 0; tstOff < dim; ++tstOff)
			{
				std::size_t m = gTst + tstOff*numUnk;
				for(unsigned basOff = 0; basOff < dim; ++basOff)
				{
					std::size_t n = gBas + basOff*numUnk;
					if(!(fixedDOFs[m] || fixedDOFs[n]))
						res += elemMat(dim*iTst + tstOff, dim*iBas + basOff)*(*pVVec)(bfRemapVec[m])*(*pVVec)(bfRemapVec[n]);
				}
			}
		}
	}
	return res;
}

std::vector<double> FEMProblem::gradCompliance(const Topologies::TOMesh* const inMesh, 
		const std::vector<std::map<std::size_t, double>>& dTOR) const
{
	std::size_t numNodes = probMesh->getNumUnknowns();
	std::size_t nelems = probMesh->getNumElements();
	std::vector<double> res(dTOR.size(), 0.);
	// Precompute uku values
	std::vector<double> ukuVec(nelems);
	for(std::size_t ielem = 0; ielem < nelems; ++ielem)
  {
		Eigen::MatrixXd elemMat = probMesh->getElementMatrix(ielem);
		ukuVec[ielem] = -elementCompliance(ielem, elemMat, numNodes);
		ukuVec[ielem] /= inMesh->getOptVal(ielem);
	}
	// Multiply with dTOR
	for(std::size_t k = 0; k < dTOR.size(); ++k)
	{
		const std::map<std::size_t, double>& curRow = dTOR[k];
		for(auto columnIt = curRow.begin(); columnIt != curRow.end(); ++columnIt)
			res[k] += ukuVec[columnIt->first]*(columnIt->second);
	}
	return res;
}

