/*
 * OdeFun.cpp
 *
 * This is a parent class which allows for different functions
 * to be used when calculating an ODE solution using RkFehl.cpp.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan
 */

#include "OdeFun.h"
#include <cmath>

using namespace std;
using namespace arma;

/**
 * Default constructor.
 */
OdeFun::OdeFun() {
	vector<double> start(1);
	start[0]=1;
	params=start;
	dim =2;
}

/**
 * Constructor with parameters provided.
 *
 * @param params - vector containing the parameters to set.
 */
OdeFun::OdeFun(vector<double> params){
	this->params=params;
	dim=2;
}

/**
 * Method to be overloaded by sub-classed for use in
 * calculating the value at the current time and state.
 *
 * @param t - current time.
 * @param x - current space state.
 * @return - output of the function.
 */
vec OdeFun::evalF(double t,vec x){
	return x;
}

OdeFun::~OdeFun() {
	// TODO Auto-generated destructor stub
}
