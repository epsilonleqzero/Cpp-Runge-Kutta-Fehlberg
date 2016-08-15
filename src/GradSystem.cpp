/*
 * GradSystem.cpp
 *
 * Implements a gradient system.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan.
 */
#include <cmath>
#include "GradSystem.h"


using namespace std;
using namespace arma;

/**
 * Default constructor with default parameters.
 */
GradSystem::GradSystem() {
	dim=2;
}

/**
 * Constructor with parameters provided.
 *
 * @param params - vector containing the parameters to set.
 */
GradSystem::GradSystem(vector<double> params){
	dim=2;
}

/**
 * Overloaded method from parent class for use in calculating
 * the value at the current time and state.
 *
 * @param t - current time.
 * @param x - current space state.
 * @return - output of the function.
 */
vec GradSystem::evalF(double t,vec x){
	vec u=zeros<vec>(2);
	u(0)=-4*x(0)*(x(0)-1)*(x(0)-0.5);
	u(1)=-2*x(1);
	return u;
}

GradSystem::~GradSystem() {

}

