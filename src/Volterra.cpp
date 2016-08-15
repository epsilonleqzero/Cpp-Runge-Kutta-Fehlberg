/*
 * Volterra.cpp
 *
 * Implements the Lotka-Volterra system.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan.
 */
#include <cmath>
#include "Volterra.h"


using namespace std;
using namespace arma;

/**
 * Default constructor with default parameters.
 */
Volterra::Volterra() {
	a=0.5;
	b=0.5;
	c=0.5;
	d=0.5;
	dim=2;
}

/**
 * Constructor with parameters provided.
 *
 * @param params - vector containing the parameters to set.
 */
Volterra::Volterra(vector<double> params){
	dim=2;
	a=params[0];
	b=params[1];
	c=params[2];
	d=params[3];
}

/**
 * Overloaded method from parent class for use in calculating
 * the value at the current time and state.
 *
 * @param t - current time.
 * @param x - current space state.
 * @return - output of the function.
 */
vec Volterra::evalF(double t,vec x){
	vec u=zeros<vec>(2);
	u(0)=a*x(0)-b*(x(1)*x(0));
	u(1)=(c*x(0)*x(1))-d*x(1);
	return u;
}

Volterra::~Volterra() {

}

