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
#include "linquad.h"
#include "UTIL/genericmaterial.h"
#include "point2d.h"
#include "point3d.h"

template <typename T> double LinearQuadrilateral<T>::xiNodes[4] = {-1., 1., 1., -1.};
template <typename T> double LinearQuadrilateral<T>::etaNodes[4] = {-1., -1., 1., 1.};

template <typename T>
LinearQuadrilateral<T>::LinearQuadrilateral(const std::vector<T*>& ptVec, const Topologies::GenericMaterial& inMat, bool debugPrint):
	Element<T>(ptVec, inMat)
{
	this->myPT = ptLinQuad;
	assert(ptVec.size() == 4);
	// Ensure correct ordering, a quad is not a simplex!
	// To determine order, find largest angle between all vectors emanating from first point
	// The vectors with the largest angle belong to opposite corners
	std::vector<T*>& locNodeVec = this->nodeVec;
	std::vector<double> angleVec(locNodeVec.size() - 1);
	for(unsigned k = 0; k < angleVec.size(); ++k)
	{
		unsigned kp1 = k + 1, kp2 = (k + 1) % angleVec.size() + 1;
		T v1 = *locNodeVec[kp1] - *locNodeVec[0], v2 = *locNodeVec[kp2] - *locNodeVec[0];
		angleVec[k] = acos(v1*v2/abs(v1)/abs(v2));
	}
	std::vector<double>::const_iterator maxIt = std::max_element(angleVec.begin(), angleVec.end());
	unsigned kmax = maxIt - angleVec.begin();
	unsigned kp1 = kmax + 1, kp2 = (kmax + 1) % angleVec.size() + 1, kp3 = (kmax + 2) % angleVec.size() + 1;
	// Recopy the node pointersa
	locNodeVec[1] = ptVec[kp2];
	locNodeVec[2] = ptVec[kp3];
	locNodeVec[3] = ptVec[kp1];
	// Sign of Jacobian is not fixed here, it will be during mesh processing if necessary
	if(debugPrint)
	{
		std::cout << "Angles: " << angleVec[0] << " " << angleVec[1] << " " << angleVec[2] << std::endl;
		std::cout << "ks: " << kp1 << " " << kp2 << " " << kp3 << std::endl;
		std::cout << "nodes: " << *locNodeVec[0] << "\n" << *locNodeVec[1] << "\n" << *locNodeVec[2] << "\n" << *locNodeVec[3] << "\n" << std::endl;
		std::cout << "nodes orig: " << *ptVec[0] << "\n" << *ptVec[1] << "\n" << *ptVec[2] << "\n" << *ptVec[3] << "\n" << std::endl;
	}
//	int bs; std::cin >> bs;
}

template <typename T>
LinearQuadrilateral<T>::~LinearQuadrilateral()
{
}

template <typename T>
T LinearQuadrilateral<T>::eval(double xi, double eta) const
{
	double x = 0., y = 0., z = 0.;
	for(unsigned k = 0; k < 4; ++k)
	{
		x += (*this->nodeVec[k])[0]*bf(k, xi, eta);
		y += (*this->nodeVec[k])[1]*bf(k, xi, eta);
		z += (*this->nodeVec[k])[2]*bf(k, xi, eta);
	}
	return T(x, y, z);
}

template <typename T>
double LinearQuadrilateral<T>::getJacobian(double xi, double eta) const
{
	T rxi = getXiVector(xi, eta), reta = getEtaVector(xi, eta);
	return rxi.crossProductNorm(reta); // Sign for 2d is handled in this function
}

template <typename T>
double LinearQuadrilateral<T>::jacdxidx(double xi, double eta) const
{
	return dydeta(xi, eta);
}

template <typename T>
double LinearQuadrilateral<T>::jacdxidy(double xi, double eta) const
{
	return -dxdeta(xi, eta);
}

template <typename T>
double LinearQuadrilateral<T>::jacdetadx(double xi, double eta) const
{
	return -dydxi(xi, eta);
}

template <typename T>
double LinearQuadrilateral<T>::jacdetady(double xi, double eta) const
{
	return dxdxi(xi, eta);
}

template <typename T>
double LinearQuadrilateral<T>::dxids(unsigned k) const
{
	if(k == 0)
		return 1;
	else if(k == 2)
		return -1;
	return 0.;
}

template <typename T>
double LinearQuadrilateral<T>::detads(unsigned k) const
{
	if(k == 1)
		return 1;
	else if(k == 3)
		return -1;
	return 0.;
}

template <typename T>
double LinearQuadrilateral<T>::dxdxi(double xi, double eta) const
{
	double res = 0.;
	for(unsigned k = 0; k < 4; ++k)
		res += this->nodeVec[k]->x*dbfdxi(k, xi, eta);
	return res;
}

template <typename T>
double LinearQuadrilateral<T>::dxdeta(double xi, double eta) const
{
	double res = 0.;
	for(unsigned k = 0; k < 4; ++k)
		res += this->nodeVec[k]->x*dbfdeta(k, xi, eta);
	return res;
}

template <typename T>
double LinearQuadrilateral<T>::dydxi(double xi, double eta) const
{
	double res = 0.;
	for(unsigned k = 0; k < 4; ++k)
		res += this->nodeVec[k]->y*dbfdxi(k, xi, eta);
	return res;
}

template <typename T>
double LinearQuadrilateral<T>::dydeta(double xi, double eta) const
{
	double res = 0.;
	for(unsigned k = 0; k < 4; ++k)
		res += this->nodeVec[k]->y*dbfdeta(k, xi, eta);
	return res;
}

template <typename T>
double LinearQuadrilateral<T>::dzdxi(double xi, double eta) const
{
	double res = 0.;
	for(unsigned k = 0; k < 4; ++k)
		res += (*this->nodeVec[k])[2]*dbfdxi(k, xi, eta); // z-component
	return res;
}

template <typename T>
double LinearQuadrilateral<T>::dzdeta(double xi, double eta) const
{
	double res = 0.;
	for(unsigned k = 0; k < 4; ++k)
		res += (*this->nodeVec[k])[2]*dbfdeta(k, xi, eta); // z-component
	return res;
}

template <typename T>
double LinearQuadrilateral<T>::jacdbfdx(unsigned k, double xi, double eta) const
{
	return dbfdxi(k, xi, eta)*jacdxidx(xi, eta) + dbfdeta(k, xi, eta)*jacdetadx(xi, eta);
}

template <typename T>
double LinearQuadrilateral<T>::jacdbfdy(unsigned k, double xi, double eta) const
{
	return dbfdxi(k, xi, eta)*jacdxidy(xi, eta) + dbfdeta(k, xi, eta)*jacdetady(xi, eta);
}

template <typename T>
T LinearQuadrilateral<T>::getXiVector(double xi, double eta) const
{
	return T(dxdxi(xi, eta), dydxi(xi, eta), dzdxi(xi, eta));
}

template <typename T>
T LinearQuadrilateral<T>::getEtaVector(double xi, double eta) const
{
	return T(dxdeta(xi, eta), dydeta(xi, eta), dzdeta(xi, eta));
}

template <typename T>
double LinearQuadrilateral<T>::getArea() const
{
	/*
	double da = 0.;
	const std::vector<T*>& locNodeVec = this->nodeVec;
	for(unsigned k = 0; k < 4; ++k)
	{
		unsigned kp1 = (k + 1) % 4;
		da += locNodeVec[k]->x*locNodeVec[kp1]->y - locNodeVec[k]->y*locNodeVec[kp1]->x;
	}
	return da/2.;*/
	// Numerically integrate the Jacobian
	double gpxi[4] = {-sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3)};
	double gpeta[4] = {-sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3)};
	double w = 1., area = 0.;
	for(unsigned kg = 0; kg < 4; ++kg)
	{
		double xi = gpxi[kg], eta = gpeta[kg], w = 1.;
		double fact = getJacobian(xi, eta);
		area += w*getJacobian(xi, eta);
	}
	return area;
}

template <typename T>
void LinearQuadrilateral<T>::switchXiEta()
{
	Element<T>::switchXiEta();
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getElemMat(const ElasticProblemType inEPT) const
{
	std::pair<double,double> cijs = Element<T>::getCijs(inEPT);
	return getDelBFElemMat(cijs.first, cijs.second);
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getDelBFElemMat(const double cii, const double cij) const
{
	// Stiffness matrix, integral of [B]^T [E] [B]
	Eigen::MatrixXd elemMat(8,8);
	elemMat.setZero();
	double mu = this->material.getParameter(Element<T>::cMu);
	// Numerical integration of stiffness matrix, 2 point integration, product rule
	double gpxi[4] = {-sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3)};
	double gpeta[4] = {-sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3)};
	// Set stiffness mat
	for(unsigned kg = 0; kg < 4; ++kg)
	{
		double xi = gpxi[kg], eta = gpeta[kg], w = 1.;
		double fact = w/getJacobian(xi, eta);
		for(int kr = 0; kr < 4; kr++)
		{
			for(int kc = 0; kc <= kr; kc++)
			{
				double elem = cii*jacdbfdx(kr, xi, eta)*jacdbfdx(kc, xi, eta) + mu*jacdbfdy(kr, xi, eta)*jacdbfdy(kc, xi, eta);
				elem *= fact;
				elemMat(kr, kc) += elem;
				if(kr != kc)
					elemMat(kc, kr) += elem;
			}
			for(int kc = 0; kc <= kr; kc++)
			{
				double elem = cij*jacdbfdx(kr, xi, eta)*jacdbfdy(kc, xi, eta) + mu*jacdbfdy(kr, xi, eta)*jacdbfdx(kc, xi, eta);
				elem *= fact;
				elemMat(kr, kc + 4) += elem;
				if(kr != kc)
					elemMat(kc + 4, kr) += elem;
			}
		}
		for(int kr = 0; kr < 4; kr++)
		{
			for(int kc = 0; kc <= kr; kc++)
			{
				double elem = cij*jacdbfdy(kr, xi, eta)*jacdbfdx(kc, xi, eta) + mu*jacdbfdx(kr, xi, eta)*jacdbfdy(kc, xi, eta);
				elem *= fact;
				elemMat(kr + 4, kc) += elem;
				if(kr != kc)
					elemMat(kc, kr + 4) += elem;
			}

			for(int kc = 0; kc <= kr; kc++)
			{
				double elem = cii*jacdbfdy(kr, xi, eta)*jacdbfdy(kc, xi, eta) + mu*jacdbfdx(kr, xi, eta)*jacdbfdx(kc, xi, eta);
				elem *= fact;
				elemMat(kr + 4, kc + 4) += elem;
				if(kr != kc)
					elemMat(kc + 4, kr + 4) += elem;
			}
		}
	}
	return elemMat;
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getMassMat() const
{
	// Integral of bf * tf
	double rho = this->material.getParameter(Element<T>::cRho);
	Eigen::MatrixXd massMat = Eigen::MatrixXd(8,8);
	for(unsigned kr = 0; kr < 8; ++kr)
		for(unsigned kc = 0; kc < 8; ++kc)
			massMat(kr,kc) = 0.;
	// Using numerical integration, will be exact for this order of basis function
	double gpxi[4] = {-sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3)};
	double gpeta[4] = {-sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3)};
	double w = 1.;
	for(unsigned kg = 0; kg < 4; ++kg)
	{
		double xi = gpxi[kg], eta = gpeta[kg];
		double jac = getJacobian(xi, eta);
		for(unsigned kr = 0; kr < 4; ++kr)
		{
			for(unsigned kc = 0; kc <= kr; ++kc)
			{
				double val = rho*w*jac*bf(kr, xi, eta)*bf(kc, xi, eta);
				massMat(kr,kc) += val;
				massMat(kr+4,kc+4) += val;
				if(kr != kc)
				{
					massMat(kc,kr) += val;
					massMat(kc+4,kr+4) += val;
				}
			}
		}
	}
	return massMat;
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getTFInt() const 
{
	// Integral of tf
	Eigen::MatrixXd fMat = Eigen::MatrixXd(8,8);
	for(int kr = 0; kr < 4; ++kr)
	{
		for(unsigned kc = 0; kc < 4; ++kc)
		{
			fMat(kr, kc + 4) = 0.;
			fMat(kr + 4, kc) = 0.;
		}
	}
	// Numerical integration, 1 point rule, exact with this bf
	double w = 4., jac = getJacobian(0., 0.);
	for(unsigned kr = 0; kr < 4; ++kr)
	{
		for(unsigned kc = 0; kc <= kr; ++kc)
		{
			double entry = w*jac*bf(kr, 0., 0.);
			fMat(kr,kc) = entry;
			fMat(kc,kr) = entry;
			fMat(kr + 4,kc + 4) = entry;
			fMat(kc + 4,kr + 4) = entry;
		}
	}
	return fMat;
/*
	std::cout << "fmmat = [";
  for(unsigned kr = 1; kr <= 8; ++kr)
  {
    for(unsigned kc = 1; kc <= 8; ++kc)
      std::cout << fMat(kr, kc) << " ";
    std::cout << std::endl;
  }
  std::cout << "]" << std::endl;
	int bs; std::cin >> bs;*/
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getLaplacianElemMat(double matVal) const
{
	// Stiffness matrix, integral of matVal * [B]^T * [B]
	Eigen::MatrixXd elemMat(4,4);
	// Numerical integration of stiffness matrix, 2 point integration, product rule
	double gpxi[4] = {-sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3)};
	double gpeta[4] = {-sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3)};
	double w = 1.;
	// Initialize elemMat
	elemMat.setZero();
	// Set stiffnes mat
	for(unsigned kg = 0; kg < 4; ++kg)
	{
		double xi = gpxi[kg], eta = gpeta[kg], w = 1.;
		double fact = matVal*w/getJacobian(xi, eta);
		for(int kr = 0; kr < 4; kr++)
		{
			for(int kc = 0; kc <= kr; kc++)
			{
				double elem = jacdbfdx(kr, xi, eta)*jacdbfdx(kc, xi, eta) + jacdbfdy(kr, xi, eta)*jacdbfdy(kc, xi, eta);
				elem *= fact;
				elemMat(kr, kc) += elem;
				if(kr != kc)
					elemMat(kc, kr) += elem;
			}
		}
	}
	return elemMat;
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getLaplacianElemMat(const Eigen::MatrixXd& matVal) const
{
	// Stiffness matrix, integral of matVal * [B]^T * [B]
	Eigen::MatrixXd elemMat(4,4);
	// Numerical integration of stiffness matrix, 2 point integration, product rule
	double gpxi[4] = {-sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3)};
	double gpeta[4] = {-sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3)};
	double w = 1.;
	// Initialize elemMat
	elemMat.setZero();
	// Set stiffnes mat
	for(unsigned kg = 0; kg < 4; ++kg)
	{
		double xi = gpxi[kg], eta = gpeta[kg], w = 1.;
		double fact = w/getJacobian(xi, eta);
		for(int kr = 0; kr < 4; kr++)
		{
			for(int kc = 0; kc <= kr; kc++)
			{
				double bfkrdx = jacdbfdx(kr, xi, eta), bfkcdx = jacdbfdx(kc, xi, eta);
				double bfkrdy = jacdbfdy(kr, xi, eta), bfkcdy = jacdbfdy(kc, xi, eta);
				double elem = matVal(0,0)*bfkrdx*bfkcdx + matVal(0,1)*bfkrdx*bfkcdy + matVal(1,0)*bfkrdy*bfkcdx + matVal(1,1)*bfkrdy*bfkcdy;
				elem *= fact;
				elemMat(kr, kc) += elem;
				if(kr != kc)
					elemMat(kc, kr) += elem;
			}
		}
	}
	return elemMat;
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getThermalExpansionMat(const ElasticProblemType inEPT, double alphaTE) const
{
	std::pair<double,double> cijs = Element<T>::getCijs(inEPT);
	double alphaIsoC = -alphaTE*(cijs.first + cijs.second);
	// Thermal expansion matrix: [B]^T * [C] * [alpha] * [N]
	Eigen::MatrixXd elemMat(8,4);
	// Numerical integration of stiffness matrix, 2 point integration, product rule
	double gpxi[4] = {-sqrt(1./3), sqrt(1./3), -sqrt(1./3), sqrt(1./3)};
	double gpeta[4] = {-sqrt(1./3), -sqrt(1./3), sqrt(1./3), sqrt(1./3)};
	double w = 1.;
	double fact = w*alphaIsoC;
	// Initialize elemMat
	elemMat.setZero();
	// Set stiffnes mat
	for(unsigned kg = 0; kg < 4; ++kg)
	{
		double xi = gpxi[kg], eta = gpeta[kg];
		for(int kr = 0; kr < 4; kr++)
			for(int kc = 0; kc < 4; kc++)
				elemMat(kr,kc) += fact*jacdbfdx(kr,xi,eta)*bf(kc,xi,eta);
		for(int kr = 0; kr < 4; kr++)
			for(int kc = 0; kc < 4; kc++)
				elemMat(kr+4,kc) += fact*jacdbfdy(kr,xi,eta)*bf(kc,xi,eta);
	}
	return elemMat;
}

template <typename T>
Eigen::VectorXd LinearQuadrilateral<T>::getInterpolationVec(double xi, double eta) const
{
	Eigen::VectorXd outVec(4);
	for(unsigned k = 0; k < 4; ++k)
		outVec(k) = bf(k, xi, eta);
	return outVec;
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getGradBF() const
{
	Eigen::MatrixXd res(2,4);
	double xi = 0., eta = 0.;
	for(int kc = 0; kc < 4; kc++)
	{
		double dbfdx = jacdbfdx(kc, xi, eta), dbfdy = jacdbfdy(kc, xi, eta);
		double jac = getJacobian(xi, eta);
		res(0,kc) = dbfdx/jac;
		res(1,kc) = dbfdy/jac;
	}
	return res;
}

template <typename T>
void LinearQuadrilateral<T>::getXiEtaLineInt(unsigned k, double& xi, double& eta) const
{
	// Integration nodes for line integrals around element contour
	if(k == 0)
	{
		xi  = 0.;
		eta = -1.;
	}
	else if(k == 1)
	{
		xi = 1.;
		eta = 0.;
	}
	else if(k == 2)
	{
		xi = 0;
		eta = 1;
	}
	else 
	{
		xi = -1;
		eta = 0.;
	}
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getLaplacianBF1() const
{
	// Line integral to compute div(grad(phi))
	// Uses weak form to give average over element (i.e. assumed constant phi over element)
	Eigen::MatrixXd outMat(1,4);
	outMat.setZero();
	// For quadrilateral, labeling is:
	//   __L3_
	//  |     |
	//L4|     | L2
	//  |_____|
	//     L1
	double area = getArea();
	for(unsigned kl = 0; kl < 4; ++kl)
	{
		double xi, eta, w = 2.;
		getXiEtaLineInt(kl, xi, eta);
		double dxds = dxdxi(xi, eta)*dxids(kl) + dxdeta(xi, eta)*detads(kl), 
			   dyds = dydxi(xi, eta)*dxids(kl) + dydeta(xi, eta)*detads(kl);
		double ds = sqrt(dxds*dxds + dyds*dyds);
		double jac = getJacobian(xi, eta);
		double nx = dyds/ds, ny = -dxds/ds;
		for(unsigned k = 0; k < 4; ++k)
			outMat(0,k) += (nx*jacdbfdx(k, xi, eta) + ny*jacdbfdy(k, xi, eta))*ds*w/jac/area;
	}
	return outMat;
}

template <typename T>
Eigen::MatrixXd LinearQuadrilateral<T>::getLaplacianBF2() const
{
	// Line integral to compute grad(grad(phi))
	// Uses weak form to give average over element (i.e. assumed constant phi over element)
	Eigen::MatrixXd outMat(4,4);
	outMat.setZero();
	// For quadrilateral, labeling is:
	//   __L3_
	//  |     |
	//L4|     | L2
	//  |_____|
	//     L1
	double area = getArea();
	for(unsigned kl = 0; kl < 4; ++kl)
	{
		double xi, eta, w = 2.;
		getXiEtaLineInt(kl, xi, eta);
		double dxds = dxdxi(xi, eta)*dxids(kl) + dxdeta(xi, eta)*detads(kl), 
			   dyds = dydxi(xi, eta)*dxids(kl) + dydeta(xi, eta)*detads(kl);
		double ds = sqrt(dxds*dxds + dyds*dyds);
		double jac = getJacobian(xi, eta);
		double nx = dyds/ds, ny = -dxds/ds;
		for(unsigned k = 0; k < 4; ++k)
		{

			outMat(0,k) += nx*jacdbfdx(k, xi, eta)*ds*w/jac/area;
			outMat(1,k) += ny*jacdbfdx(k, xi, eta)*ds*w/jac/area;
			outMat(2,k) += nx*jacdbfdy(k, xi, eta)*ds*w/jac/area;
			outMat(3,k) += ny*jacdbfdy(k, xi, eta)*ds*w/jac/area;
		}
	}
	return outMat;
}

template class LinearQuadrilateral<Point2D>;
template class LinearQuadrilateral<Point3D>;
