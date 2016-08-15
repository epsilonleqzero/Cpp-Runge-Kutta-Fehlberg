/*
 * RkFehl.cpp
 *
 * This function runs the Runge-Kutta Fehlberg method on an
 * ODE selected by the string given as an input. It uses
 * adaptive time-stepping to create the output matrix.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan
 */

#include "RkFehl.h"

using namespace std;
using namespace arma;

RkFehl::RkFehl() {
	h = 0.5;
	dim = 1;
	vector<double> pars(2);
	pars[0] = 1.0;
	pars[1] = 1.0;
}

/**
 * Constructor to run the RK45 method. Takes the inputs then calculates the solution
 * to the ODE.
 *
 * @param params - Vector containing parameters for the ODE functions.
 * @param fun - String containing which function to run.
 * @param tf - Double representing the final time.
 * @param h - Normal step-size value. Used for error.
 * @param x - Initial conditions.
 */
RkFehl::RkFehl(vector<double> params, string fun, double tf, double h, vec x) {
	// Choose the correct ODE to solve.
	if ((fun.find("Duff") != string::npos)
			|| (fun.find("duff") != string::npos)) {
		ode = new Duffing(params);
		dim = ode->dim;
	} else if ((fun.find("Van") != string::npos)
			|| (fun.find("van") != string::npos)) {
		ode = new VanDerPol(params);
		dim = ode->dim;
	} else if ((fun.find("Lor") != string::npos)
			|| (fun.find("lor") != string::npos)) {
		ode = new Lorenz(params);
		dim = ode->dim;
	} else if ((fun.find("Ross") != string::npos)
			|| (fun.find("ross") != string::npos)) {
		ode = new Rossler(params);
		dim = ode->dim;
	} else if ((fun.find("Grad") != string::npos)
			|| (fun.find("grad") != string::npos)) {
		ode = new GradSystem(params);
		dim = ode->dim;
	} else if ((fun.find("clinic") != string::npos)
			|| (fun.find("Clinic") != string::npos)) {
		ode = new Heteroclinic(params);
		dim = ode->dim;
	} else if ((fun.find("olterra") != string::npos)
			|| (fun.find("Lotka") != string::npos)) {
		ode = new Volterra(params);
		dim = ode->dim;
	} else {
		ode = new GradSystem(params);
		dim = ode->dim;
	}
	// Store parameters from input.
	this->h = h;
	vec y1 = x;
	t.push_back(0);
	y.push_back(y1);
	double tcurr = t[0];
	double hc = h;
	int i = 1;
	// Setup error checking and truncation error for adaptive
	// time stepping.
	double epsloc = pow(h, 4);
	double kappa = 0.84;
	while (tcurr < tf) {
		// Calculate the test step.
		vector<vec> yv = testStep(tcurr, y[i - 1], hc);
		// Evaluate the error at the current step.
		double epsi = (norm(yv[0] - yv[1]));
		int j = 0;
		// Set optimal h.
		if (epsi < epsloc) {
			hc = hc * kappa * pow(epsloc / epsi, 0.25);
		}
		// Re-run if the error is too large on the optimal time step
		// value.
		while (epsi > epsloc && j < 30) {
			j++;
			hc = hc * kappa * pow(epsloc / epsi, 0.25);
			yv = testStep(tcurr, y[i - 1], hc);
			epsi = (norm(yv[0] - yv[1]));
		}
		// Save solution at current step.
		tcurr = t[i - 1] + hc;
		t.push_back(tcurr);
		y.push_back(yv[1]);
		i++;
	}
	// Erase memory no longer needed.
	delete ode;
}


/**
 * Calculate the next step for the RK45 method. Unused for now.
 *
 * @param t - current time.
 * @param x - current values.
 * @return - solution at next time.
 */
vec RkFehl::nextStep(double t, vec x) {
	double hc = h;
	vector<vec> yv = testStep(t, x, hc);
	return x;
}

/**
 * Calculates the RK45 method on the current step. This
 * test is used to multiple times to ensure that the error is
 * kept low.
 *
 * @param t - current time.
 * @param x - current values.
 * @param hc - current time-step.
 * @return both of the function values to test for error.
 */
vector<vec> RkFehl::testStep(double t, vec x, double hc) {
	// RK4 vector.
	vec y4 = zeros<vec>(dim);
	// RK5 vector.
	vec y5 = zeros<vec>(dim);

	// Calculate values according to the RK45 algorithm.
	vec k1 = ode->evalF(t, x) * hc;
	vec k2 = ode->evalF(t + hc / 4.0, x + k1 / 4.0) * hc;
	vec k3 = ode->evalF(t + (3.0 * hc / 8.0),
			x + (3.0 * k1 / 32.0) + (9.0 * k2 / 32.0)) * hc;
	vec k4 = ode->evalF(t + (12.0 * hc / 13.0),
			x + (1932.0 * k1 / 2197.0) - (7200.0 * k2 / 2197.0)
					+ (7296.0 * k3 / 2197.0)) * hc;
	vec k5 = ode->evalF(t + hc,
			x + (439.0 * k1 / 216.0) - (8.0) * k2 + (3680.0 * k3 / 513.0)
					- (845.0 * k4 / 4104.0)) * hc;
	vec k6 = ode->evalF(t + hc / 2,
			x - (8.0 / 27.0) * k1 + (2 * k2) - (3544.0 * k3 / 2565.0)
					+ (1859.0 * k4 / 4104.0) - (11.0 / 40.0) * k5) * hc;
	// Calculate the next step values.
	y4 = x + (25 * k1 / 216.0) + (1408.0 * k3 / 2565.0) + (2197.0 * k4 / 4101.0)
			- (k5 / 5);
	y5 = x + (16 * k1 / 135.0) + (6656.0 * k3 / 12825.0)
			+ (28561.0 * k4 / 56430.0) - (9.0 * k5 / 50.0) + (2.0 * k6 / 55.0);
	// Setup output.
	vector<vec> ret(2);
	ret[0] = y4;
	ret[1] = y5;
	return ret;
}

RkFehl::~RkFehl() {

}

