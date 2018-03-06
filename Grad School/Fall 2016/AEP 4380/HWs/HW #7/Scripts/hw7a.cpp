/* AEP 4380 HW#7a
  Time-dependent Schrodinger equation
  Run on core i7 with g++ 5.4.0 (Ubuntu)
  Gianfranco Grillo 11/01/2016
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>
#define ARRAYT_BOUNDS_CHECK
#include "arrayt.hpp"
#include <complex>
typedef complex<double> CMPLX;

using namespace std;

int main() {
    double deltax = 0.1, maxx = 1000.0, x;
    int j, xlength;
    xlength = maxx/deltax;
    arrayt<double> xlist(xlength);
    double pot(double);
    CMPLX initwf(double), wf0;
    ofstream fp;
    fp.open( "task1.dat" );

    for (j = 0; j < xlength; j++) {
        x = j*deltax;
        wf0 = initwf(x);
        fp << x << setw(15) << wf0.real() << setw(15) << wf0.imag() << setw(15) << norm(wf0) << setw(15) << pot(x) << endl;
    }
    fp.close();
    return( EXIT_SUCCESS );
}

CMPLX initwf(double x) {
    double L = 1000.0, s = 20.0, re, im, y, m;
    y = exp(-(x-0.3*L)*(x-0.3*L)/(s*s));
    re = cos(x)*y;
    im = sin(x)*y;
    return CMPLX(re, im);
}

double pot(double x) {
    double L = 1000.0, omv = 7.0, V0 = 4.05;
    return (V0/(1.0+exp((0.5*L-x)/omv)));
}
