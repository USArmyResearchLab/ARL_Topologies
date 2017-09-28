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

#include "trilinhex.h"

double TriLinHex::xi1Nodes[8] = {-1.,  1.,  1., -1., -1.,  1., 1., -1.};
double TriLinHex::xi2Nodes[8] = {-1., -1.,  1.,  1., -1., -1., 1.,  1.};
double TriLinHex::xi3Nodes[8] = {-1., -1., -1., -1.,  1.,  1., 1.,  1.};
double TriLinHex::xi1Gauss1[8] = {-sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3)};
double TriLinHex::xi2Gauss1[8] = {-sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3), -sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3)};
double TriLinHex::xi3Gauss1[8] = {-sqrt(1./3), -sqrt(1./3), -sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3), sqrt(1./3), sqrt(1./3)};
double TriLinHex::wGauss1[8] = {1., 1., 1., 1., 1., 1., 1., 1.};

TriLinHex::TriLinHex(const std::vector<Point3D*>& inPtVec, const Topologies::GenericMaterial& inMaterial):
	Hexahedron(inPtVec, inMaterial)
{
	itsCellType = ctTriLinHex;
}

void TriLinHex::switchXi2Xi3()
{
	Hexahedron::switchXi2Xi3();
}

void TriLinHex::justL(double xi1, double xi2, double xi3, Point3D& l1, Point3D& l2, Point3D& l3) const
{
	Eigen::MatrixXd rmat = dxdxi(xi1, xi2, xi3);
	Point3D rXi1(rmat(0,0), rmat(0,1), rmat(0,2)), rXi2(rmat(1,0), rmat(1,1), rmat(1,2)), rXi3(rmat(2,0), rmat(2,1), rmat(2,2));
	l1 = rXi3;
	l2 = -rXi2;
	l3 = rXi1;
}

void TriLinHex::jacL(double xi1, double xi2, double xi3, double& jac, 
	Point3D& l1, Point3D& l2, Point3D& l3) const
{
	Eigen::MatrixXd rmat = dxdxi(xi1, xi2, xi3);
	Point3D rXi1(rmat(0,0), rmat(0,1), rmat(0,2)), rXi2(rmat(1,0), rmat(1,1), rmat(1,2)), rXi3(rmat(2,0), rmat(2,1), rmat(2,2));
	l1 = rXi3;		
	l2 = -rXi2;
	l3 = rXi1;
	// Jacobian given by Eq. (42) in GWP paper and translates to:
	Point3D temp = rXi2.crossProduct(rXi3);
	jac = rXi1*temp;
}

Point3D TriLinHex::pos(double xi1, double xi2, double xi3) const
{
	Point3D tmp(0., 0., 0.);
	for(unsigned k = 0; k < nodeVec.size(); ++k)
		tmp += bf(k, xi1, xi2, xi3)*(*nodeVec[k]);
	return tmp;
}

Eigen::MatrixXd TriLinHex::getElemMat() const
{
	Eigen::MatrixXd E = elasticityTensor();
	Eigen::MatrixXd B = bmat(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
	double dJ = getJacobian(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
	Eigen::MatrixXd tmp =  wGauss1[0]*B.transpose()*E*B*dJ;
	for(unsigned k = 1; k < nGauss1; ++k)
	{
		B = bmat(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		dJ = getJacobian(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		tmp += wGauss1[k]*B.transpose()*E*B*dJ;
	}
	return tmp;
}

Eigen::MatrixXd TriLinHex::getMassMat() const
{
	double rho = itsMaterial.getParameter(0);
	Eigen::MatrixXd massMat(24,24);
	massMat.setZero();
	for(unsigned kg = 0; kg < nGauss1; ++kg)
	{
		double jac = getJacobian(xi1Gauss1[kg], xi2Gauss1[kg], xi3Gauss1[kg]);
		for(unsigned kr = 0; kr < 8; ++kr)
		{
			unsigned rowOff = kr*3;
			for(unsigned kc = 0; kc < 8; ++kc)
			{
				unsigned colOff = kc*3;
				double val = rho*wGauss1[kg]*jac*bf(kr, xi1Gauss1[kg], xi2Gauss1[kg], xi3Gauss1[kg])*bf(kc, xi1Gauss1[kg], xi2Gauss1[kg], xi3Gauss1[kg]);
				massMat(rowOff + 0, colOff + 0) += val;
				massMat(rowOff + 1, colOff + 1) += val;
				massMat(rowOff + 2, colOff + 2) += val;
			}
		}
	}
	return massMat;
}

Eigen::MatrixXd TriLinHex::getTFInt() const
{
	// Integral of tf
	Eigen::MatrixXd tfMat(24,24);
	tfMat.setZero();
	for(unsigned kg = 0; kg < nGauss1; ++kg)
	{
		double jac = getJacobian(xi1Gauss1[kg], xi2Gauss1[kg], xi3Gauss1[kg]);
		for(unsigned kr = 0; kr < 8; ++kr)
		{
			unsigned rowOff = kr*3;
			for(unsigned kc = 0; kc < 8; ++kc)
			{
				unsigned colOff = kc*3;
				double val = wGauss1[kg]*jac*bf(kr, xi1Gauss1[kg], xi2Gauss1[kg], xi3Gauss1[kg]);
				tfMat(rowOff + 0, colOff + 0) += val;
				tfMat(rowOff + 1, colOff + 1) += val;
				tfMat(rowOff + 2, colOff + 2) += val;
			}
		}
	}
	return tfMat;
}

double TriLinHex::volumeIntegral() const
{
	// Integrate Jacobian
	double vol = 0.;
	for(unsigned k = 0; k < nGauss1; ++k)
		vol += wGauss1[k]*getJacobian(xi1Gauss1[k],xi2Gauss1[k],xi3Gauss1[k]);
	return vol;
}

double TriLinHex::getJacobian(double xi1, double xi2, double xi3) const
{
	Eigen::MatrixXd rmat = dxdxi(xi1, xi2, xi3);
	Point3D rXi1(rmat(0,0), rmat(0,1), rmat(0,2)), rXi2(rmat(1,0), rmat(1,1), rmat(1,2)), rXi3(rmat(2,0), rmat(2,1), rmat(2,2));
	Point3D temp = rXi2.crossProduct(rXi3);
	return rXi1*temp;
}

Eigen::MatrixXd TriLinHex::dxdxi(double xi1, double xi2, double xi3) const
{
	Eigen::MatrixXd rmat(3,8);
	for(unsigned k = 0; k < 8; ++k)
	{
		rmat(0,k) = nodeVec[k]->x;
		rmat(1,k) = nodeVec[k]->y;
		rmat(2,k) = nodeVec[k]->z;
	}
	return rmat*dNdxi(xi1, xi2, xi3);
}

Eigen::MatrixXd TriLinHex::elasticityTensor() const
{
	double lambda = itsMaterial.getParameter(1), mu = itsMaterial.getParameter(2);
	double cii = lambda + 2.*mu, cij = lambda;
	Eigen::MatrixXd cijkl(6,6);
	cijkl.setZero();
	cijkl(0,0) = cii;
	cijkl(1,1) = cii;
	cijkl(2,2) = cii;
	cijkl(0,1) = cij;
	cijkl(0,2) = cij;
	cijkl(1,0) = cij;
	cijkl(1,2) = cij;
	cijkl(2,0) = cij;
	cijkl(2,1) = cij;
	cijkl(3,3) = mu;
	cijkl(4,4) = mu;
	cijkl(5,5) = mu;
	return cijkl;
}

Eigen::MatrixXd TriLinHex::bmat(double xi1, double xi2, double xi3) const
{
	// Gives the derivative of the basis functions, in a form convenient for computing strain/stress
	Eigen::MatrixXd dndxmat = dNdx(xi1, xi2, xi3);
	Eigen::MatrixXd b(6,24);
	b.setZero();
	for(unsigned k = 0; k < 8; ++k)
	{
		unsigned koff = 3*k;
		b(0, 0 + koff) = dndxmat(k, 0);
		b(1, 1 + koff) = dndxmat(k, 1);
		b(2, 2 + koff) = dndxmat(k, 2);
		b(3, 0 + koff) = dndxmat(k, 1);
		b(3, 1 + koff) = dndxmat(k, 0);
		b(4, 1 + koff) = dndxmat(k, 2);
		b(4, 2 + koff) = dndxmat(k, 1);
		b(5, 0 + koff) = dndxmat(k, 2);
		b(5, 2 + koff) = dndxmat(k, 0);
	}
	return b;
}

Eigen::MatrixXd TriLinHex::dNdx(double xi1, double xi2, double xi3) const
{
	Eigen::MatrixXd dndximat = dNdxi(xi1, xi2, xi3);
	Eigen::MatrixXd invJ = inverseJac(xi1, xi2, xi3);
	return dndximat*invJ;
}

Eigen::MatrixXd TriLinHex::dNdxi(double xi1, double xi2, double xi3) const
{
	Eigen::MatrixXd gradN(8,3);
	for(unsigned k = 0; k < 8; ++k)
	{
		gradN(k,0) = dbfdxi1(k, xi1, xi2, xi3);
		gradN(k,1) = dbfdxi2(k, xi1, xi2, xi3);
		gradN(k,2) = dbfdxi3(k, xi1, xi2, xi3);
	}
	return gradN;
}

Eigen::MatrixXd TriLinHex::inverseJac(double xi1, double xi2, double xi3) const
{
	Eigen::MatrixXd invJ(3,3);
	double dJ = getJacobian(xi1, xi2, xi3);
	Eigen::MatrixXd J = dxdxi(xi1, xi2, xi3);
	invJ(0,0) = (J(1,1)*J(2,2) - J(1,2)*J(2,1))/dJ;
	invJ(0,1) = (J(0,2)*J(2,1) - J(0,1)*J(2,2))/dJ;
	invJ(0,2) = (J(0,1)*J(1,2) - J(0,2)*J(1,1))/dJ;

	invJ(1,0) = (J(1,2)*J(2,0) - J(1,0)*J(2,2))/dJ;
	invJ(1,1) = (J(0,0)*J(2,2) - J(0,2)*J(2,0))/dJ;
	invJ(1,2) = (J(0,2)*J(1,0) - J(0,0)*J(1,2))/dJ;

	invJ(2,0) = (J(1,0)*J(2,1) - J(1,1)*J(2,0))/dJ;
	invJ(2,1) = (J(0,1)*J(2,0) - J(0,0)*J(2,1))/dJ;
	invJ(2,2) = (J(0,0)*J(1,1) - J(0,1)*J(1,0))/dJ;

	return invJ;
}

Eigen::MatrixXd TriLinHex::getLaplacianElemMat(double matVal) const
{
	Eigen::MatrixXd B = dNdx(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
	double dJ = getJacobian(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
	Eigen::MatrixXd tmp =  wGauss1[0]*B*matVal*B.transpose()*dJ;
	for(unsigned k = 1; k < nGauss1; ++k)
	{
		B = dNdx(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		dJ = getJacobian(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		tmp += wGauss1[k]*B*matVal*B.transpose()*dJ;
	}
	return tmp;
}

Eigen::MatrixXd TriLinHex::getLaplacianElemMat(const Eigen::MatrixXd& matVal) const
{
	Eigen::MatrixXd B = dNdx(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
	double dJ = getJacobian(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
	Eigen::MatrixXd tmp =  wGauss1[0]*B*matVal*B.transpose()*dJ;
	for(unsigned k = 1; k < nGauss1; ++k)
	{
		B = dNdx(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		dJ = getJacobian(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		tmp += wGauss1[k]*B*matVal*B.transpose()*dJ;
	}
	return tmp;
}

Eigen::MatrixXd TriLinHex::getThermalExpansionMat(double alphaTE) const
{
	Eigen::MatrixXd B = bmat(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
  Eigen::VectorXd N = getInterpolationVec(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
  Eigen::MatrixXd C = elasticityTensor();
  // Set up thermal expansion tensor and include Jacobian
  Eigen::VectorXd alphaMat(6);
  alphaMat.setZero();
  alphaMat(0) = -alphaTE;
  alphaMat(1) = -alphaTE;
  alphaMat(2) = -alphaTE;
	// Integrate
	double dJ = getJacobian(xi1Gauss1[0], xi2Gauss1[0], xi3Gauss1[0]);
	Eigen::MatrixXd tmp = (wGauss1[0]*dJ)*B.transpose()*C*alphaMat*N.transpose();
	for(unsigned k = 1; k < nGauss1; ++k)
	{
		B = bmat(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		N = getInterpolationVec(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		dJ = getJacobian(xi1Gauss1[k], xi2Gauss1[k], xi3Gauss1[k]);
		tmp += (wGauss1[k]*dJ)*B.transpose()*C*alphaMat*N.transpose();
	}
  return tmp;
}

Eigen::VectorXd TriLinHex::getInterpolationVec(double xi, double eta, double zeta) const
{
	Eigen::VectorXd outVec(8);
	for(unsigned k = 0; k < 8; ++k)
		outVec(k) = bf(k, xi, eta, zeta);
	return outVec;
}

Eigen::MatrixXd TriLinHex::getGradBF() const
{
	return Eigen::MatrixXd();
}

Eigen::MatrixXd TriLinHex::getLaplacianBF1() const
{
	return Eigen::MatrixXd();
}

Eigen::MatrixXd TriLinHex::getLaplacianBF2() const
{
	return Eigen::MatrixXd();
}
