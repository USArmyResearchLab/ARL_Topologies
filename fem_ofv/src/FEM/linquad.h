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

#ifndef LINQUAD_H
#define LINQUAD_H

#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "element.h"

namespace Topologies{
class GenericMaterial;
}

//! LinearQuadrilateral implements an Element as a flat edged quadrilateral with bi-linear basis functions
/*! This class implements bi-linear basis functions, meaning that they are not strictly linear.  They are 
 *  complete to first order, but, assuming that the area coordinates are \xi and \eta, they contain the 
 *  higher order term \xi \eta.  
 *  As a quadrilateral is not a simplex, it requires some preprocessing to ensure that the nodes are in
 *  a standard order.  This is done in the constructor so the nodes passed to the constructor can be
 *  in any order.
 */
template <class T>
class LinearQuadrilateral : public Element<T>
{
public:
	//! Constructor taking a vector of defining nodes and material properties.
	/*! The nodes can be in any order, the constructor will reorder them */
	LinearQuadrilateral(const std::vector<T*>& ptVec, const Topologies::GenericMaterial& inMat, bool debugPrint = false);
	virtual ~LinearQuadrilateral();

	//! @name Element matrix computation functions
	//@{
	virtual T eval(double xi, double eta) const;
	virtual Eigen::MatrixXd getElemMat(const ElasticProblemType inEPT) const;
	virtual Eigen::MatrixXd getLaplacianElemMat(double matVal) const;
	virtual Eigen::MatrixXd getLaplacianElemMat(const Eigen::MatrixXd& matVal) const;
	virtual Eigen::MatrixXd getMassMat() const;
	virtual Eigen::MatrixXd getTFInt() const;
	virtual Eigen::MatrixXd getThermalExpansionMat(const ElasticProblemType inEPT, double alphaTE) const;
	virtual Eigen::VectorXd getInterpolationVec(double xi, double eta) const;
	virtual Eigen::MatrixXd getGradBF() const;
	virtual Eigen::MatrixXd getLaplacianBF1() const;
	virtual Eigen::MatrixXd getLaplacianBF2() const;
	//@}
	//! @name Geometric functions
	//@{
	virtual double getArea() const;
	virtual double getJacobian(double xi, double eta) const;
	virtual T getXiVector(double xi, double eta) const;
	virtual T getEtaVector(double xi, double eta) const;
	//@}
	//! @name Topological functions
	//@{
	virtual int getOrientation() const;
	virtual void switchXiEta();
	//@}
private:
	Eigen::MatrixXd getDelBFElemMat(const double cii, const double cij) const;
	double bf(unsigned k, double xi, double eta) const {return 0.25*(1 + xi*xiNodes[k])*(1 + eta*etaNodes[k]);}
	double dbfdxi(unsigned k, double xi, double eta) const {return 0.25*xiNodes[k]*(1 + eta*etaNodes[k]);}
	double dbfdeta(unsigned k, double xi, double eta) const {return 0.25*(1 + xi*xiNodes[k])*etaNodes[k];}
	double jacdxidx(double xi, double eta) const;
	double jacdxidy(double xi, double eta) const;
	double jacdetadx(double xi, double eta) const;
	double jacdetady(double xi, double eta) const;
	double dxdxi(double xi, double eta) const;
	double dxdeta(double xi, double eta) const;
	double dydxi(double xi, double eta) const;
	double dydeta(double xi, double eta) const;
	double dzdxi(double xi, double eta) const;
	double dzdeta(double xi, double eta) const;
	double jacdbfdx(unsigned k, double xi, double eta) const;
	double jacdbfdy(unsigned k, double xi, double eta) const;
	double dxids(unsigned k) const;
	double detads(unsigned k) const;
	void getXiEtaLineInt(unsigned k, double& xi, double& eta) const;

	static double xiNodes[4], etaNodes[4];
};

template <typename T> inline
int LinearQuadrilateral<T>::getOrientation() const
{
	return getArea() > 0 ? 1 : -1;
}

#endif
