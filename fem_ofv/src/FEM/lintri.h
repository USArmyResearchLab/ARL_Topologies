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

#ifndef LINETRIANGLE_H
#define LINETRIANGLE_H

#include <Eigen/Dense>
#include "UTIL/topologiesdefs.h"
#include "element.h"

namespace Topologies{
class GenericMaterial;
}

//! LinearTriangle implements an Element as a flat edged triangle with linear basis functions
/*! Since a triangle is a 2d simplex, this is the simplest implementation of Element.  
*/
template <class T>
class LinearTriangle : public Element<T>
{
public:
	//! Constructor taking a vector of 3 defining nodes and material properties
	LinearTriangle(const std::vector<T*>& ptVec, const Topologies::GenericMaterial& inMat);
	virtual ~LinearTriangle();

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
	void getDiffJacBFs(Eigen::VectorXd& dNdxJac, Eigen::VectorXd& dNdyJac) const;

	T r00, rxi, reta;
	double da;
};

template <typename T> 
double LinearTriangle<T>::getArea() const
{
	return fabs(da)/2.;
}

template <typename T>
int LinearTriangle<T>::getOrientation() const
{
	return da > 0 ? 1 : -1;
}

template <typename T>
double LinearTriangle<T>::getJacobian(double xi, double eta) const
{
	return da;
}

template <typename T>
T LinearTriangle<T>::getXiVector(double xi, double eta) const
{
	return rxi;
}

template <typename T>
T LinearTriangle<T>::getEtaVector(double xi, double eta) const
{
	return reta;
}

#endif
