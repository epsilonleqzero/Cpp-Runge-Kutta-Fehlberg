
/*
 * Volterra.h
 *
 *  Created on: Aug 6, 2016
 *      Author: Devils
 */
#include<vector>
#include<armadillo>
#include "OdeFun.h"

#ifndef VOLTERRA_H_
#define VOLTERRA_H_

class Volterra: public OdeFun {
public:
	Volterra();
	virtual ~Volterra();
	Volterra(std::vector<double> params);
	arma::vec evalF(double t,arma::vec x);
	int dim;
private:
	std::vector<double> params;
	double a;
	double b;
	double c;
	double d;
	};
#endif /* VOLTERRA_H_ */
