/*
 * Rossler.h
 *
 *  Created on: Aug 6, 2016
 *      Author: Devils
 */
#include<vector>
#include<armadillo>
#include "OdeFun.h"

#ifndef ROSSLER_H_
#define ROSSLER_H_

class Rossler: public OdeFun {
public:
	Rossler();
	virtual ~Rossler();
	Rossler(std::vector<double> params);
	arma::vec evalF(double t,arma::vec x);
	int dim;
private:
	std::vector<double> params;
	double a;
	double b;
	double c;
	};
#endif /* DUFFING_H_ */
