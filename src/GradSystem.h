/*
 * GradSystem.h
 *
 *  Created on: Aug 6, 2016
 *      Author: Devils
 */
#include<vector>
#include<armadillo>

#ifndef GRADSYSTEM_H_
#define GRADSYSTEM_H_

#include "OdeFun.h"

class GradSystem: public OdeFun {
public:
	GradSystem();
	virtual ~GradSystem();
	GradSystem(std::vector<double> params);
	arma::vec evalF(double t,arma::vec x);
	int dim;
private:
	std::vector<double> params;
	};
#endif /* GRADSYSTEM_H_ */
