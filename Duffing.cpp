/*
 * Duffing.cpp
 *
 * Implements the Duffing Oscillator.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan.
 */
#include <cmath>
#include "Duffing.h"


using namespace std;
using namespace arma;

/**
 * Default constructor with default parameters.
 */
Duffing::Duffing() {
	delta=0.5;
	gamma=0.5;
	w=0.5;
    bta=0.5;
	dim=2;
}

/**
 * Constructor with parameters provided.
 *
 * @param params - vector containing the parameters to set.
 */
Duffing::Duffing(vector<double> params){
	dim=2;
    bta=params[0];
	gamma=params[2];
	w=params[3];
	delta=params[1];
}

/**
 * Overloaded method from parent class for use in calculating
 * the value at the current time and state.
 *
 * @param t - current time.
 * @param x - current space state.
 * @return - output of the function.
 */
vec Duffing::evalF(double t,vec x){
	vec u=zeros<vec>(2);
	u(0)=x(1);
	u(1)=-bta*x(0)-pow(x(0),3.0)-delta*x(1)+gamma*cos(w*t);
	return u;
}

Duffing::~Duffing() {

}

