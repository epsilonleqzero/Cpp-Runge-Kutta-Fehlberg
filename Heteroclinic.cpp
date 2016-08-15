/*
 * Heteroclinic.cpp
 *
 * Implements a heteroclinic cycle.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan.
 */
#include <cmath>
#include "Heteroclinic.h"


using namespace std;
using namespace arma;

/**
 * Default constructor with default parameters.
 */
Heteroclinic::Heteroclinic() {
	dim=2;
	epsilon=0.0;
	mu1=0.1;
	mu2=0.1;
}

/**
 * Constructor with parameters provided.
 *
 * @param params - vector containing the parameters to set.
 */
Heteroclinic::Heteroclinic(vector<double> params){
	dim=2;
	epsilon=params[0];
	mu1=params[1];
	mu2=params[2];
}

/**
 * Overloaded method from parent class for use in calculating
 * the value at the current time and state.
 *
 * @param t - current time.
 * @param x - current space state.
 * @return - output of the function.
 */
vec Heteroclinic::evalF(double t,vec x){
	vec u=zeros<vec>(2);
	u(0)=x(1)+epsilon*(mu1*x(0)+mu2*x(1));
	u(1)=-x(0)+pow(x(0),3.0);
	return u;
}

Heteroclinic::~Heteroclinic() {

}

