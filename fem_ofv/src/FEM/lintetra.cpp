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

#include "lintetra.h"

LinTetra::LinTetra(const std::vector<Point3D*>& inPtVec, const Topologies::GenericMaterial& inMaterial):
Tetrahedron(inPtVec, inMaterial),
	r(*inPtVec[0]),
	rXi1(*inPtVec[1] - *inPtVec[0]),
	rXi2(*inPtVec[2] - *inPtVec[0]),
	rXi3(*inPtVec[3] - *inPtVec[0])
{
	itsCellType = ctLinTetra;
}

void LinTetra::switchXi2Xi3()
{
	std::swap(rXi2, rXi3);
	Tetrahedron::switchXi2Xi3();
}

void LinTetra::justL(double xi1, double xi2, double xi3, Point3D& l1, Point3D& l2, Point3D& l3) const
{
	assert(xi1 >= 0. && xi2 >= 0. && xi3 >= 0. && xi1 + xi2 + xi3 <= 1.);
	l1 = rXi3;		// l1 is L12 directed component (GWP paper)
	l2 = -rXi2;		// l2 is L13 directed component (GWP paper)
	l3 = rXi1;		// l3 is L23 directed component (GWP paper)
}

void LinTetra::jacL(double xi1, double xi2, double xi3, double& jac, 
	Point3D& l1, Point3D& l2, Point3D& l3) const
{
	Point3D temp;
	assert(xi1 >= 0. && xi2 >= 0. && xi3 >= 0. && xi1 + xi2 + xi3 <= 1.);
	l1 = rXi3;		
	l2 = -rXi2;		
	l3 = rXi1;

	// Jacobian given by Eq. (42) in GWP paper and translates to:
	temp = -l2.crossProduct(l1);
	jac = l3*temp;
}

Point3D LinTetra::pos(double xi1, double xi2, double xi3) const
{
	assert(xi1 >= 0. && xi2 >= 0. && xi3 >= 0. && xi1 + xi2 + xi3 <= 1.);
	return r + xi1*rXi1 + xi2*rXi2 + xi3*rXi3;
}

Eigen::MatrixXd LinTetra::getElemMat() const
{
	Eigen::MatrixXd B = bmat(), E = elasticityTensor();
	double dJ = getJacobian();
	Eigen::MatrixXd tmp = B.transpose()*E*B*(dJ/6.);
	return tmp;
}

Eigen::MatrixXd LinTetra::getMassMat() const
{
	double rho = itsMaterial.getParameter(0);
	double entry = rho*getJacobian()/120.;
	Eigen::MatrixXd massMat(12,12);
	massMat.setZero();
	for(unsigned kr = 0; kr < 4; ++kr)
	{
		unsigned rowOff = kr*3;
		for(unsigned kc = 0; kc < 4; ++kc)
		{
			unsigned colOff = kc*3;
			massMat(rowOff + 0, colOff + 0) = entry;
			massMat(rowOff + 1, colOff + 1) = entry;
			massMat(rowOff + 2, colOff + 2) = entry;
		}
	}
	for(unsigned kr = 0; kr < 12; kr++)
		massMat(kr,kr) *= 2.;
	return massMat;
}

Eigen::MatrixXd LinTetra::getTFInt() const
{
	// Integral of tf
	double entry = getJacobian()/24.;
	Eigen::MatrixXd tfMat(12,12);
	tfMat.setZero();
	for(int kr = 0; kr < 4; ++kr)
	{
		unsigned rowOff = kr*3;
		for(unsigned kc = 0; kc < 4; ++kc)
		{
			unsigned colOff = kc*3;
			tfMat(rowOff + 0, colOff + 0) = entry;
			tfMat(rowOff + 1, colOff + 1) = entry;
			tfMat(rowOff + 2, colOff + 2) = entry;
		}
	}
	return tfMat;
}

double LinTetra::volumeIntegral() const
{
	return getJacobian()/6.;
}

double LinTetra::getJacobian() const
{
	Point3D temp = rXi2.crossProduct(rXi3);
	return rXi1*temp;
}

Eigen::MatrixXd LinTetra::dxdxi() const
{
	Eigen::MatrixXd dxmat(3,3);
	dxmat(0,0) = rXi1.x;
	dxmat(0,1) = rXi2.x;
	dxmat(0,2) = rXi3.x;
	dxmat(1,0) = rXi1.y;
	dxmat(1,1) = rXi2.y;
	dxmat(1,2) = rXi3.y;
	dxmat(2,0) = rXi1.z;
	dxmat(2,1) = rXi2.z;
	dxmat(2,2) = rXi3.z;
	return dxmat;
}

Eigen::MatrixXd LinTetra::elasticityTensor() const
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

Eigen::MatrixXd LinTetra::bmat() const
{
	// Gives the derivative of the basis functions, in a form convenient for computing strain/stress
	Eigen::MatrixXd dndxmat = dNdx();
	Eigen::MatrixXd b(6,12);
	b.setZero();
	for(unsigned k = 0; k < 4; ++k)
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

Eigen::MatrixXd LinTetra::dNdx() const
{
	Eigen::MatrixXd dndximat = dNdxi();
	Eigen::MatrixXd invJ = inverseJac();
	Eigen::MatrixXd tmp = dndximat*invJ;
	return tmp;
}

Eigen::MatrixXd LinTetra::dNdxi() const
{
	Eigen::MatrixXd gradN(4,3);
	// BFs: b1 = 1 - xi - eta - zeta
	//      b2 = xi
	//      b3 = eta
	//      b4 = zeta
	gradN(0,0) = -1.;
	gradN(0,1) = -1.;
	gradN(0,2) = -1.;
	gradN(1,0) =  1.;
	gradN(1,1) =  0.;
	gradN(1,2) =  0.;
	gradN(2,0) =  0.;
	gradN(2,1) =  1.;
	gradN(2,2) =  0.;
	gradN(3,0) =  0.;
	gradN(3,1) =  0.;
	gradN(3,2) =  1.;
	return gradN;
}

Eigen::MatrixXd LinTetra::inverseJac() const
{
	Eigen::MatrixXd invJ(3,3);
	double dJ = getJacobian();
	Eigen::MatrixXd J = dxdxi();
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

Eigen::MatrixXd LinTetra::getLaplacianElemMat(double matVal) const
{
	Eigen::MatrixXd B = dNdx();
	double dJ = getJacobian();
	Eigen::MatrixXd tmp =  B*B.transpose()*(dJ*matVal/6.);
	return tmp;
}

Eigen::MatrixXd LinTetra::getLaplacianElemMat(const Eigen::MatrixXd& matVal) const
{
	Eigen::MatrixXd B = dNdx();
	double dJ = getJacobian();
	Eigen::MatrixXd tmp =  B*matVal*B.transpose()*(dJ/6.);
	return tmp;
}

Eigen::MatrixXd LinTetra::getThermalExpansionMat(double alphaTE) const
{
	Eigen::MatrixXd B = bmat();
	Eigen::VectorXd N = getInterpolationVec(0.25,0.25,0.25);
	Eigen::MatrixXd C = elasticityTensor();
	// Set up thermal expansion tensor and include Jacobian
	Eigen::VectorXd alphaMat(6);
	alphaMat.setZero();
	double dJ = getJacobian()/6.;
	alphaMat(0) = -alphaTE*dJ;
	alphaMat(1) = -alphaTE*dJ;
	alphaMat(2) = -alphaTE*dJ;
	return B.transpose()*C*alphaMat*N.transpose();
}

Eigen::VectorXd LinTetra::getInterpolationVec(double xi, double eta, double zeta) const
{
	Eigen::VectorXd outVec(4);
	outVec(0) = 1. - xi - eta - zeta;
	outVec(1) = xi;
	outVec(2) = eta;
	outVec(3) = zeta;
	return outVec;
}

Eigen::MatrixXd LinTetra::getGradBF() const
{
	return Eigen::MatrixXd();
}

Eigen::MatrixXd LinTetra::getLaplacianBF1() const
{
	return Eigen::MatrixXd();
}

Eigen::MatrixXd LinTetra::getLaplacianBF2() const
{
	return Eigen::MatrixXd();
}
