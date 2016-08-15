/*
 * Lorenz.cpp
 *
 * Implements the Lorenz equations.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan.
 */
#include <cmath>
#include "Lorenz.h"

using namespace std;
using namespace arma;

/**
 * Default constructor with default parameters.
 */
Lorenz::Lorenz() {
	sigma=1.0;
	rho=1.0;
	beta=1.0;
	dim=3;
}

/**
 * Constructor with parameters provided.
 *
 * @param params - vector containing the parameters to set.
 */
Lorenz::Lorenz(vector<double> params){
	dim=3;
	beta=params[0];
	rho=params[1];
	sigma=params[2];
}

 /**
  * Overloaded method from parent class for use in calculating
  * the value at the current time and state.
  *
  * @param t - current time.
  * @param x - current space state.
  * @return - output of the function.
  */
vec Lorenz::evalF(double t,vec x){
	vec u=zeros<vec>(dim);
	u(0)=sigma*(x(1)-x(0));
	u(1)=x(0)*(rho-x(2))-x(1);
	u(2)=x(0)*x(1)-beta*x(2);
	return u;
}

Lorenz::~Lorenz() {
	// TODO Auto-generated destructor stub
}

