/* AEP 4380 Final Project
  Solving the Stellar Structure Equations
  Run on core i5 with g++ 5.4.0
  Gianfranco Grillo 12/03/2016
*/

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>
#include "arrayt.hpp"
using namespace std;

int main() {
    int npts = 1e5, i, j, n = 4, nsteps = 1000, minpos;
    double K = 3.841e14, pi = 4.0 * atan(1.0), M = 1.9891e33, h, hold, start, stop, stepw, m = 0.0, centralp;
    bool cond = true;
    arrayt<double> initvar(n), ans(n), err(nsteps), finalvar(n), presslist(nsteps);
    void autosteprk4(arrayt<double>&, arrayt<double>&, int, double&, double, void(arrayt<double>&, double, arrayt<double>&, double&));
    void stellar(arrayt<double>&, double, arrayt<double>&, double&);
    int minloc(arrayt<double>&);
    ofstream fp;
    start = 5e16;
    stop = 1e18;
    while (cond == true) {
        stepw = (stop-start)/nsteps;
        h = M/npts;
        for (i = 0; i < nsteps; i++) {
            initvar(0) = 1;
            presslist(i) = start + i*stepw;
            initvar(1) = presslist(i);
            initvar(2) = 2e7;
            initvar(3) = 1;
            m = 0.0;
            while (m < M) {
                hold = h;
                autosteprk4(initvar, finalvar, n, h, m, stellar);
                if (finalvar(1) < 0.0 || finalvar(2) < 0.0) break;
                if (h > hold || h == hold) { // means error is low enough, result is correct
                    for (j = 0; j < n; j++) { initvar(j) = finalvar(j); }
                    m = m + hold;
                }
            }
            err(i) = abs((m-M)/M)+abs(initvar(1)/(start+i*stepw)+abs(initvar(2)/1e7));
            cout << err(i) << endl;
        }
        minpos = minloc(err);
        if (minloc == 0) { // no minimum found, increase initial guess range
            start = start/4;
            stop = 4*stop;
        }
        else { cond = false; }
    }
    initvar(0) = 1;
    initvar(1) = presslist(minpos);
    initvar(2) = 1e7;
    initvar(3) = 1;
    fp.open("solution.dat");
    m = 0.0;
    while (m < M) {
        hold = h;
        autosteprk4(initvar, finalvar, n, h, m, stellar);
        if (finalvar(1) < 0.0 || finalvar(2) < 0.0) break;
        if (h > hold || h == hold) { // means error is low enough, result is correct
            for (i = 0; i < n; i++) { initvar(i) = finalvar(i); }
            fp << m << setw(15) << initvar(0) << setw(15) << initvar(1) << setw(15) << initvar(2) << setw(15) << initvar(3) << endl;
            m = m + hold;
        }
    }
    fp.close();
    cout << "Mass = " << m << setw(15) << "Radius = " << initvar(0) << "       " << "Central pressure = " << presslist(minpos) << "    " << "Luminosity = " << initvar(3) << endl;
    return (EXIT_SUCCESS);
}

void autosteprk4(arrayt<double>& initvar, arrayt<double>& highfinalvar, int n, double& h,
				 double t0, void functions(arrayt<double>&, double, arrayt<double>&, double&)) {
	int i;
	arrayt<double> scale(n), error(n), lowfinalvar(n), k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n), scalek(n), temp(n);
    double maxerror = 1e-8, deltamax;
	double maximum(arrayt<double>&);
	double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0/9.0, c6 = 1.0, c7 = 1.0;
	double a21 = 0.2, a31 = 3.0/40.0, a32 = 9.0/40.0, a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0;
	double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0;
	double a61 = 9017.0/3168.0, a62 = -355.0/33.0, a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0;
	double b1 = 35.0/384.0, b3 = 500.0/1113.0, b4 = 125.0/192.0, b5 = -2187.0/6784.0, b6 = 11.0/84.0;
	double d1 = 5179.0/57600.0, d3 = 7571.0/16695.0, d4 = 393.0/640.0, d5 = -92097.0/339200.0, d6 = 187.0/2100.0, d7 = 1.0/40.0;
	functions(initvar, t0, k1, h); // generate k1
	for (i = 0; i < n; i++) {
		temp(i) = initvar(i) + a21*k1(i);
	}
	functions(temp, t0 + c2*h, k2, h);
	for (i = 0; i < n; i++) {
		temp(i) = initvar(i) + a31*k1(i) + a32*k2(i);
	}
	functions(temp, t0 + c3*h, k3, h);
	for (i = 0; i < n; i++) {
		temp(i) = initvar(i) + a41*k1(i) + a42*k2(i) + a43*k3(i);
	}
	functions(temp, t0 + c4*h, k4, h);
	for (i = 0; i < n; i++) {
		temp(i) = initvar(i) + a51*k1(i) + a52*k2(i) + a53*k3(i) + a54*k4(i);
	}
	functions(temp, t0 + c5*h, k5, h);
	for (i = 0; i < n; i++) {
		temp(i) = initvar(i) + a61*k1(i) + a62*k2(i) + a63*k3(i) + a64*k4(i) + a65*k5(i);
	}
	functions(temp, t0 + c6*h, k6, h);
	for (i = 0; i < n; i++) {
		highfinalvar(i) = initvar(i) + b1*k1(i) + b3*k3(i) + b4*k4(i) + b5*k5(i) + b6*k6(i);
	}
	functions(highfinalvar, t0 + h, k7, h);
	for (i = 0; i < n; i++) {
		lowfinalvar(i) = initvar(i) + d1*k1(i) + d3*k3(i) + d4*k4(i) + d5*k5(i) + d6*k6(i) + d7*k7(i);
	}
	functions(initvar, t0 + h, scalek, h); // generate scale
	for (i = 0; i < n; i++) {
		scale(i) = abs(initvar(i)) + abs(scalek(i)) + 0.01;
	}
	for (i = 0; i < n; i++) {
		error(i) = abs((highfinalvar(i) - lowfinalvar(i))/scale(i));
	}
	deltamax = maximum(error);
	if (5*deltamax < maxerror) {
		h = h*pow(maxerror/deltamax, 0.2); // if error is too small, increase h
		return;
	}
	else if (deltamax > maxerror && abs(h) > 1e-12) { // if error is too big, decrease h
		h = h/3.0;
		return;
	}
	else {
		return;
	}
}

void stellar(arrayt<double>& y, double m, arrayt<double>& k, double& h) {
    double K = 3.841e14, pi = 4.0 * atan(1.0), r, P, G = 6.6732e-8, mp = 1.672e-24, nu, eps0, kb = 1.380e-16, mu = 0.6, T;
    r = y(0);
    P = y(1);
    T = y(2);
    k(0) = h*pow(K,3.0/4.0)/(4*pi*r*r*pow(P, 3.0/4.0));
    k(1) = -h*G*m/(4*pi*r*r*r*r);
    k(2) = -h*G*m*pow(K,0.75)*mu*mp/(16*pi*r*r*r*r*kb*pow(P,0.75));
    eps0 = -3.90476190e-23*T*T*T+1.29666667e-14*T*T-2.68666667e-07*T+1.49904762;
    nu = -1.209*log(T)+23.98;
    k(3) = h*eps0*P*mu*mp*pow(T,nu-1)/kb;
	return;
}

double maximum(arrayt<double>& array) {
	double max;
	int i, n = array.n();
	max = array(0);
	for (i = 0; i < n; i++) {
		if (array(i) > max) {
			max = array(i);
		}
	}
	return max;
}

int minloc(arrayt<double>& array) {
    double m;
	int i, n = array.n(), minpos;
	m = array(0);
    minpos = 0;
	for (i = 0; i < n; i++) {
		if (array(i) < m && isnan(array(i))==false) {
			m = array(i);
            minpos = i;
		}
	}
    return minpos;
}
