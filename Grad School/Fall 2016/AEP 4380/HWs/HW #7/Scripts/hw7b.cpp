/* AEP 4380 HW#7b
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
#include "arrayt.hpp"
#include <complex>
typedef complex<double> CMPLX;
typedef arrayt<CMPLX> arrayc;

using namespace std;

int main() {
    double deltax = 0.1, deltat, maxt = 5e-14, maxx = 1000.0, x, t, h = 6.5821e-16, m = 5.699e-32, w, par1, par2;
    int i, j, N, tlength, counter=0;

    cout << "Please enter sampling time in seconds: "; // ask for sampling time input
    cin >> deltat; 

    N = maxx/deltax;
    w = 2.0*deltax*deltax/deltat;
    par1 = 2.0*m*w/h;
    par2 = 2.0*m*deltax*deltax/(h*h);
    tlength = maxt/deltat+1;
    arrayt<double> xlist(N), tlist(tlength);
    arrayc alist(N-1), blist(N), clist(N-1), dlist(N), psilist(N), clistor(N-1);
    void dlistcreator(arrayc&, arrayc&, double, double, arrayt<double>&);
    void tridiag(arrayc&, arrayc&, arrayc&, arrayc&, arrayc&);
    double pot(double);
    CMPLX initwf(double), wf0;
    ofstream fp;
    fp.open("task2.dat");

    for (j = 0; j < N; j++) { xlist(j) = j*deltax; }

    for (j = 0; j < tlength; j++) { tlist(j) = j*deltat; }

    for (j = 0; j < N; j++) { // initialize arrays
        x = xlist(j);
        wf0 = initwf(x);
        psilist(j) = wf0;
        blist(j) = CMPLX((-2.0 - par2*pot(x)), par1);
        if (j != N-1) { alist(j) = CMPLX(1.0, 0.0); clist(j) = CMPLX(1.0, 0.0); }
    }
    for (j = 0; j < tlength; j++) {
        t = tlist(j);
        if (counter == (tlength-1)/5 || counter == 2*(tlength-1)/5 || counter == 3*(tlength-1)/5 || counter == 4*(tlength-1)/5 || counter == (tlength-1)) {
            for (i = 0; i < N; i++) {
                fp << t << setw(15) << xlist(i) << setw(15) << norm(psilist(i)) << endl; //output to file only in five instances
            }
        }
        dlistcreator(dlist, psilist, par1, par2, xlist);
        tridiag(alist, blist, clist, dlist, psilist);
        counter++;
    }
    fp.close();
    return( EXIT_SUCCESS );
}

void dlistcreator(arrayc& dlist, arrayc& psilist, double par1, double par2, arrayt<double>& xlist) {
    double x;    
    int N = psilist.n(), j=0;
    double pot(double);
    for (j = 0; j < N; j++) {
        x = xlist(j);
        if (j == 0) {
            dlist(j) = CMPLX((2.0 + par2*pot(x)), par1)*psilist(0)-psilist(1);
        }
        else if (j == (N-1)) {
            dlist(j) = CMPLX((2.0 + par2*pot(x)), par1)*psilist(j)-psilist(j-1);
        }
        else {
            dlist(j) = CMPLX((2.0 + par2*pot(x)), par1)*psilist(j)-psilist(j-1)-psilist(j+1);
        }
    }
    return;
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

void tridiag(arrayc &alist, arrayc &blist, arrayc &clist, arrayc &dlist, arrayc &psilist) {
    int j, N;
    CMPLX bet;
    N = alist.n();
    arrayc gam(N);
    psilist(0) = dlist(0)/(bet=blist(0));
    for (j=1; j<N; j++){
        gam(j) = clist(j-1)/bet;
        bet = blist(j)-alist(j)*gam(j);
        psilist(j) = (dlist(j)-alist(j)*psilist(j-1))/bet;
    }
    for (j=(N-2); j >= 0; j--)
        psilist(j) -= gam(j+1)*psilist(j+1);
    return;
}
