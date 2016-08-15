/*
 * Rossler.cpp
 *
 * Implements the Rossler attractor.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan.
 */
#include <cmath>
#include "Rossler.h"

using namespace std;
using namespace arma;

/**
 * Default constructor with default parameters.
 */
Rossler::Rossler() {
	a=0.2;
	b=0.2;
	c=5.7;
	dim=3;
}

/**
 * Constructor with parameters provided.
 *
 * @param params - vector containing the parameters to set.
 */
Rossler::Rossler(vector<double> params){
	dim=3;
	a=params[0];
	b=params[1];
	c=params[2];
}

/**
 * Overloaded method from parent class for use in calculating
 * the value at the current time and state.
 *
 * @param t - current time.
 * @param x - current space state.
 * @return - output of the function.
 */
vec Rossler::evalF(double t,vec x){
	vec u=zeros<vec>(dim);
	u(0)=-x(1)-x(2);
	u(1)=x(0)+a*x(1);
	u(2)=b+x(2)*(x(0)-c);
	return u;
}

Rossler::~Rossler() {
	// TODO Auto-generated destructor stub
}

