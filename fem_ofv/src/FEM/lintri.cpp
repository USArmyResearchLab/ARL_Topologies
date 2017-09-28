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

#include <assert.h>
#include "lintri.h"
#include "UTIL/genericmaterial.h"
#include "point2d.h"
#include "point3d.h"

template <typename T>
LinearTriangle<T>::LinearTriangle(const std::vector<T*>& ptVec, const Topologies::GenericMaterial& inMat):
	Element<T>(ptVec, inMat)
{
	this->myPT = ptLinTri;
	assert(ptVec.size() == 3);
	r00 = *ptVec[0];
	rxi = *ptVec[1] - *ptVec[0];
	reta = *ptVec[2] - *ptVec[0];
	da = rxi.crossProductNorm(reta);
}

template <typename T>
LinearTriangle<T>::~LinearTriangle()
{
}

template <typename T>
T LinearTriangle<T>::eval(double xi, double eta) const
{
	T tmp = r00 + rxi*xi + reta*eta;
	return tmp;
}

template <typename T>
void LinearTriangle<T>::switchXiEta()
{
	std::swap(rxi, reta);
	da = rxi.x*reta.y - reta.x*rxi.y;
	Element<T>::switchXiEta();
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getElemMat(const ElasticProblemType inEPT) const
{
	std::pair<double,double> cijs = Element<T>::getCijs(inEPT);
	return getDelBFElemMat(cijs.first, cijs.second);
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getDelBFElemMat(const double cii, const double cij) const
{
	Eigen::MatrixXd elemMat(6,6);
	Eigen::VectorXd xVec(3);
	Eigen::VectorXd yVec(3);
	getDiffJacBFs(yVec, xVec);
	double fact = 2.*da;
	double mu = this->material.getParameter(Element<T>::cMu);
	for(unsigned kr = 0; kr < 3; kr++)
	{
		for(unsigned kc = 0; kc <= kr; kc++)
		{
			double elem = cii*yVec(kr)*yVec(kc) + mu*xVec(kr)*xVec(kc);
			elem /= fact;
			elemMat(kr, kc) = elem;
			elemMat(kc, kr) = elem;
		}

		for(unsigned kc = 0; kc <= kr; kc++)
		{
			double elem = cij*yVec(kr)*xVec(kc) + mu*xVec(kr)*yVec(kc);
			elem /= fact;
			elemMat(kr, kc + 3) = elem;
			elemMat(kc + 3, kr) = elem;
		}
	}
	for(unsigned kr = 0; kr < 3; kr++)
	{
		for(unsigned kc = 0; kc <= kr; kc++)
		{
			double elem = cij*xVec(kr)*yVec(kc) + mu*yVec(kr)*xVec(kc);
			elem /= fact;
			elemMat(kr + 3, kc) = elem;
			elemMat(kc, kr + 3) = elem;
		}

		for(unsigned kc = 0; kc <= kr; kc++)
		{
			double elem = cii*xVec(kr)*xVec(kc) + mu*yVec(kr)*yVec(kc);
			elem /= fact;
			elemMat(kr + 3, kc + 3) = elem;
			elemMat(kc + 3, kr + 3) = elem;
		}
	}
	return elemMat;
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getMassMat() const
{
	// Integral of bf * tf
	double rho = this->material.getParameter(Element<T>::cRho);
	double entry = rho*da/24.;
	Eigen::MatrixXd massMat = Eigen::MatrixXd(6,6);
	for(int kr = 0; kr < 3; ++kr)
		for(unsigned kc = 3; kc < 6; ++kc)
			massMat(kr,kc) = 0.;
	for(int kr = 3; kr < 6; ++kr)
		for(unsigned kc = 0; kc < 3; ++kc)
			massMat(kr,kc) = 0.;
	for(int kr = 0; kr < 3; ++kr)
		for(unsigned kc = 0; kc < 3; ++kc)
			massMat(kr,kc) = entry;
	for(int kr = 3; kr < 6; ++kr)
		for(unsigned kc = 3; kc < 6; ++kc)
			massMat(kr,kc) = entry;
	for(int kr = 0; kr < 3; kr++)
	{
		massMat(kr,kr) *= 2.;
		massMat(kr+3,kr+3) *= 2.;
	}
	return massMat;
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getTFInt() const 
{
	// Integral of tf
	double entry = da/6.;
	Eigen::MatrixXd fMat = Eigen::MatrixXd(6,6);
	for(int kr = 0; kr < 3; ++kr)
		for(unsigned kc = 3; kc < 6; ++kc)
			fMat(kr,kc) = 0.;
	for(int kr = 3; kr < 6; ++kr)
		for(unsigned kc = 0; kc < 3; ++kc)
			fMat(kr,kc) = 0.;
	for(int kr = 0; kr < 3; ++kr)
		for(unsigned kc = 0; kc < 3; ++kc)
			fMat(kr,kc) = entry;
	for(int kr = 3; kr < 6; ++kr)
		for(unsigned kc = 4; kc < 6; ++kc)
			fMat(kr,kc) = entry;
	return fMat;
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getLaplacianElemMat(double matVal) const
{
	Eigen::MatrixXd elemMat(3,3);
	Eigen::VectorXd dNdxJac(3);
	Eigen::VectorXd dNdyJac(3);
	getDiffJacBFs(dNdxJac, dNdyJac);
	double fact = matVal/(2.*da);
	for(int kr = 0; kr < 3; kr++)
	{
		for(int kc = 0; kc <= kr; kc++)
		{
			double elem = dNdxJac(kr)*dNdxJac(kc) + dNdyJac(kr)*dNdyJac(kc);
			elem *= fact;
			elemMat(kr, kc) = elem;
			elemMat(kc, kr) = elem;
		}
	}
	return elemMat;
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getLaplacianElemMat(const Eigen::MatrixXd& matVal) const
{
	Eigen::MatrixXd elemMat(3,3);
	Eigen::VectorXd dNdxJac(3);
	Eigen::VectorXd dNdyJac(3);
	getDiffJacBFs(dNdxJac, dNdyJac);
	double fact = 2.*da;
	for(int kr = 0; kr < 3; kr++)
	{
		for(int kc = 0; kc <= kr; kc++)
		{
			double elem = matVal(1,1)*dNdxJac(kr)*dNdxJac(kc) + matVal(1,2)*dNdxJac(kr)*dNdyJac(kc)
					  + matVal(2,1)*dNdyJac(kr)*dNdxJac(kc) + matVal(2,2)*dNdyJac(kr)*dNdyJac(kc);
			elem /= fact;
			elemMat(kr, kc) = elem;
			elemMat(kc, kr) = elem;
		}
	}
	return elemMat;
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getThermalExpansionMat(const ElasticProblemType inEPT, double alphaTE) const
{
	std::pair<double,double> cijs = Element<T>::getCijs(inEPT);
	double alphaIsoC = -alphaTE*(cijs.first + cijs.second);
	Eigen::MatrixXd elemMat(6,3);
	Eigen::MatrixXd dbfs = getGradBF();
	double bfint = 1./6.;
	for(unsigned kr = 0; kr < 3; ++kr)
		for(unsigned kc = 0; kc < 3; ++kc)
			elemMat(kr,kc) = da*alphaIsoC*bfint*dbfs(0,kr);
	for(unsigned kr = 0; kr < 3; ++kr)
		for(unsigned kc = 0; kc < 3; ++kc)
			elemMat(kr+3,kc) = da*alphaIsoC*bfint*dbfs(1,kr);
	return elemMat;
}

template <typename T>
Eigen::VectorXd LinearTriangle<T>::getInterpolationVec(double xi, double eta) const
{
	Eigen::VectorXd outVec(3);
	outVec(0) = 1. - xi - eta;
	outVec(1) = xi;
	outVec(2) = eta;
	return outVec;
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getGradBF() const
{
	// Returns a matrix containing the gradient of the basis functions
	// First row is x-component, second is y-component, columns are the corresponding basis functions
	Eigen::MatrixXd res(2,3);
	Eigen::VectorXd dNdxJac(3);
	Eigen::VectorXd dNdyJac(3);
	getDiffJacBFs(dNdxJac, dNdyJac);
	for(int kc = 0; kc < 3; kc++)
	{
		res(0,kc) = dNdxJac(kc)/da;
		res(1,kc) = dNdyJac(kc)/da;
	}
	return res;
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getLaplacianBF1() const
{
	// Line integral to compute div(grad(phi))
	// Uses weak form to give average over element (i.e. assumed constant phi over element)
	Eigen::MatrixXd outMat(1,3);
	Eigen::VectorXd dNdxJac(3);
	Eigen::VectorXd dNdyJac(3);
	getDiffJacBFs(dNdxJac, dNdyJac);
	// For triangle, labeling is:
	//  |\
	//L2| \ L1
	//  |__\
	//   L3
	// And integral is computed counter-clockwise
	// L1
	double dxds = reta.x - rxi.x, dyds = reta.y - rxi.y;
	double ds = sqrt(dxds*dxds + dyds*dyds);
	double nx = dyds/ds, ny = -dxds/ds;
	double fact = 1./da/getArea();
	for(unsigned k = 0; k < 3; ++k)
		outMat(0,k) = (nx*dNdxJac(k) + ny*dNdyJac(k))*ds*fact;
	// L2
	dxds = -reta.x;
	dyds = -reta.y;
	ds = sqrt(dxds*dxds + dyds*dyds);
	nx =  dyds/ds;
	ny = -dxds/ds;
	for(unsigned k = 0; k < 3; ++k)
		outMat(0,k) += (nx*dNdxJac(k) + ny*dNdyJac(k))*ds*fact;
	// L3
	dxds = rxi.x;
	dyds = rxi.y;
	ds = sqrt(dxds*dxds + dyds*dyds);
	nx =  dyds/ds;
	ny = -dxds/ds;
	for(unsigned k = 0; k < 3; ++k)
		outMat(0,k) += (nx*dNdxJac(k) + ny*dNdyJac(k))*ds*fact;
	return outMat;
}

template <typename T>
Eigen::MatrixXd LinearTriangle<T>::getLaplacianBF2() const
{
	// Line integral to compute grad(grad(phi))
	// Uses weak form to give average over element (i.e. assumed constant phi over element)
	Eigen::MatrixXd outMat(4,3);
	Eigen::VectorXd dNdxJac(3);
	Eigen::VectorXd dNdyJac(3);
	getDiffJacBFs(dNdxJac, dNdyJac);
	// For triangle, labeling is:
	//  |\
	//L2| \ L1
	//  |__\
	//   L3

	// L1
	double dxds = reta.x - rxi.x, dyds = reta.y - rxi.y;
	double ds = sqrt(dxds*dxds + dyds*dyds);
	double nx = dyds/ds, ny = -dxds/ds;
	double fact = 2.*ds/da/da;
	for(unsigned k = 0; k < 3; ++k)
	{
		outMat(0,k) = nx*dNdxJac(k)*fact;
		outMat(1,k) = ny*dNdxJac(k)*fact;
		outMat(2,k) = nx*dNdyJac(k)*fact;
		outMat(3,k) = ny*dNdyJac(k)*fact;
	}
	// L2
	dxds = -reta.x;
	dyds = -reta.y;
	ds = sqrt(dxds*dxds + dyds*dyds);
	nx =  dyds/ds;
	ny = -dxds/ds;
	fact = 2.*ds/da/da;
	for(unsigned k = 0; k < 3; ++k)
	{
		outMat(0,k) += nx*dNdxJac(k)*fact;
		outMat(1,k) += ny*dNdxJac(k)*fact;
		outMat(2,k) += nx*dNdyJac(k)*fact;
		outMat(3,k) += ny*dNdyJac(k)*fact;
	}
	// L3
	dxds = rxi.x;
	dyds = rxi.y;
	ds = sqrt(dxds*dxds + dyds*dyds);
	nx =  dyds/ds;
	ny = -dxds/ds;
	fact = 2.*ds/da/da;
	for(unsigned k = 0; k < 3; ++k)
	{
		outMat(0,k) += nx*dNdxJac(k)*fact;
		outMat(1,k) += ny*dNdxJac(k)*fact;
		outMat(2,k) += nx*dNdyJac(k)*fact;
		outMat(3,k) += ny*dNdyJac(k)*fact;
	}
	return outMat;
}

template <typename T>
void LinearTriangle<T>::getDiffJacBFs(Eigen::VectorXd& dNdxJac, Eigen::VectorXd& dNdyJac) const
{
	dNdxJac = Eigen::VectorXd(3);
	dNdyJac = Eigen::VectorXd(3);
	dNdxJac(0) = rxi.y - reta.y;
	dNdxJac(1) = reta.y;
	dNdxJac(2) = -rxi.y;
	dNdyJac(0) = reta.x - rxi.x;
	dNdyJac(1) = -reta.x;
	dNdyJac(2) = rxi.x;
}

// Class instantiations
template class LinearTriangle<Point2D>;
template class LinearTriangle<Point3D>;
