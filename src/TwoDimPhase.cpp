#include <iostream>
#include <armadillo>
#include "RkFehl.h"
#include "TwoDimPhase.h"

using namespace std;
using namespace arma;

/**
 * This is a Java JNI interface file. It runs calls the libraries present
 * and passes the parameters.
 *
 *
 * input double vector has 5 fields:
 *
 * h - index 0.
 * Parameters:
 * a - index 1.
 * b - index 2.
 * c - index 3.
 * d - index 4.
 * Final time:
 * tf - index 5.
 * Initial Conditions:
 * x0 - index 6.
 * y0 - index 7.
 * z0 - index 8.
 *
 *
 * @param *env - Java environment (provided)
 * @param jobj - Java object (provided)
 * @param jarray - Java double array to extract parameters from.
 * @param pptype - String representing the name of the ODE to solve.
 * @return returns a double array representing the solution to the ODE.
 */
JNIEXPORT jdoubleArray JNICALL Java_net_tedkwan_javafemjni_TwoDimPhase_rk4(
		JNIEnv *env, jobject jobj, jdoubleArray jarray, jstring pptype) {
	unsigned int i;
	// Extract information from Java.
	jboolean isCopy1;
	jboolean isCopy2;
	const char *s = env->GetStringUTFChars(pptype, &isCopy2);
	string pptypec(s);
	jdouble* srcArrayElems = env->GetDoubleArrayElements(jarray, &isCopy1);
	// Setup parameters needed to call the ODE solver.
	vec x;
	int r=3;
	// Place parameters into a vector.
	vector<double> pars(4);
	double h = srcArrayElems[0];
	pars[0] = srcArrayElems[1];
	pars[1] = srcArrayElems[2];
	pars[2] = srcArrayElems[3];
	pars[3] = srcArrayElems[4];
	double tf = srcArrayElems[5];
	// Chose function to run.
	if (pptypec.find("oren") != string::npos ||
			pptypec.find("ossle") != string::npos) {
		x = zeros<vec>(3);
		x(0) = srcArrayElems[6];
		x(1) = srcArrayElems[7];
		x(2) = srcArrayElems[8];
		r=4;
	} else {
		x = zeros<vec>(2);
		x(0) = srcArrayElems[6];
		x(1) = srcArrayElems[7];
		r=3;
	}
	// Run the Runge-Kutta Fehlberg method (adaptive).
	RkFehl testrk(pars, pptypec, tf, h, x);
	// Release memory of input string.
	if (isCopy2 == JNI_TRUE) {
		env->ReleaseStringUTFChars(pptype, s);
	}
	// Initialize output vector.
	vector<vec> y = testrk.y;
	vector<double> t = testrk.t;
	vector<double> res;
	unsigned int k = t.size();
	// Add data to output array.
	for (int j = 0; j < r; j++) {
		for (i = 0; i < k; i++) {
			// Check to see if we add the time, and
			// check to see the dimensions of the solution.
			if (j < (r-1)) {
				res.push_back(y[i](j));
			} else {
				res.push_back(t[i]);
			}
		}
	}

	// Release memory for input array.
	if (isCopy1 == JNI_TRUE) {
		env->ReleaseDoubleArrayElements(jarray, srcArrayElems, JNI_ABORT);
	}

	// Setup output array.
	int lenres = res.size();
	jdoubleArray result = env->NewDoubleArray(lenres);
	env->SetDoubleArrayRegion(result, 0, lenres, &res[0]);

	return result;
}
