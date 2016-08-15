/*
 * Heteroclinic.h
 *
 *  Created on: Aug 6, 2016
 *      Author: Devils
 */
#include<vector>
#include<armadillo>

#ifndef HETEROCLINIC_H_
#define HETEROCLINIC_H_

#include "OdeFun.h"

class Heteroclinic: public OdeFun {
public:
	Heteroclinic();
	virtual ~Heteroclinic();
	Heteroclinic(std::vector<double> params);
	arma::vec evalF(double t,arma::vec x);
	int dim;
private:
	std::vector<double> params;
	double epsilon;
	double mu1;
	double mu2;
	};
#endif /* HETEROCLINIC_H_ */
