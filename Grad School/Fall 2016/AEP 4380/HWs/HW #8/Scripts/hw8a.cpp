/* AEP 4380 HW#8a
  Least Squares Curve Fitting
  Run on core i7 with g++ 5.4.0 (Ubuntu)
  Gianfranco Grillo 11/07/2016
*/

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "arrayt.hpp"
#include "nr3.h"
#define ARRAYT_BOUNDS_CHECK

using namespace std;

int main() {
    int i, npts = 607, npar = 7, j, k, t;
    double pi = 4.0 * atan(1.0), Flk, bl, chisqr;
    arrayt<int> tlist(npts);
    arrayt<double> co2list(npts), blist(npar), fit(npts), sqrerrlist(npts), reslist(npts);
    arrayt<double> Fmatrix(npar, npar), funcm(npts,npar);
    void readfile(arrayt<int>&, arrayt<double>&); // create t, co2, arrays from file
    void gaussj(arrayt<double>&, arrayt<double>&);
    ofstream fp, sp;
    fp.open("output1.dat");
    sp.open("output2.dat");
    readfile(tlist, co2list);
    for (i = 0; i < npts; i++) { // initialize function arrays
        t = tlist(i);
        funcm(i,0) = sin(pi*t/6.0);
        funcm(i,1) = cos(pi*t/6.0);
        funcm(i,2) = sin(pi*t/3.0);
        funcm(i,3) = cos(pi*t/3.0);
        funcm(i,4) = t*t;
        funcm(i,5) = t; // don't really need this but makes life easier
        funcm(i,6) = 1.0;
        sqrerrlist(i) = (0.002*co2list(i))*(0.002*co2list(i)); // create error array
    }
    for (i = 0; i < npar; i++) for (j = 0; j < i+1; j++) { // only need to calculate half of matrix
        Flk = 0;
        for (k = 0; k < npts; k++) {
            Flk = Flk + funcm(k,i)*funcm(k,j)/sqrerrlist(k);
        }
        Fmatrix(i,j) = Flk;
        if (i != j) {
            Fmatrix(j,i) = Flk;
        }
    }
    for (i = 0; i < npar; i++) { // generate b vector
        bl = 0;
        for (k = 0; k < npts; k++) {
            bl = bl + co2list(k)*funcm(k,i)/sqrerrlist(k);
        }
        blist(i) = bl;
    }
    gaussj(Fmatrix, blist);
    chisqr = 0.0;
    for (i = 0; i < npts; i++) {
        fit(i) = 0;
        for (j = 0; j < npar; j++) {
            fit(i) = fit(i) + blist(j)*funcm(i, j);
            }
        reslist(i) = co2list(i)-fit(i);
        chisqr = chisqr + reslist(i)*reslist(i)/sqrerrlist(i);
        fp << tlist(i) << setw(15) << co2list(i) << setw(15) << fit(i) << setw(15) << reslist(i) << endl;
    }
    chisqr = chisqr/(npts-npar);
    cout << chisqr << endl;
    for (i = 0; i < npar; i++) {
        sp << Fmatrix(i,i) << setw(15) << blist(i) << endl;
    }
    fp.close();
    return (EXIT_SUCCESS);
}

void readfile(arrayt<int>& tlist, arrayt<double>& co2list) {
    int i, j, npts, year, t, nval;
    double co2, ymin, ymax;
    string cline;
    vector<double> x, y;
    ifstream fp;
    fp.open("maunaloa.co2.txt");

    for (i=0; i<15; i++) getline(fp,cline);
    t = 0;
    npts = 0;
    ymin = 1000.0;
    ymax = -ymin;
    for (i=0; i<70; i++) {
        fp >> year;
        for (j=0; j<12; j++) {
            fp >> co2;
            if (co2 > 0.0) {
                x.push_back(t);
                y.push_back(co2);
                if (y[npts] > ymax) ymax = y[npts];
                if (y[npts] < ymin) ymin = y[npts];
                tlist(npts) = t;
                co2list(npts) = co2;
                npts++;
            }
            t += 1;
        }
        if (year >=2008) break;
        getline(fp,cline);
    }
    return;
}

void gaussj(arrayt<double> &a, arrayt<double> &b) {
    int i, icol, irow, j, k, l, ll, n=a.n1();
    double big, dum, pivinv;
    arrayt<int> indxc(n), indxr(n), ipiv(n);
    for (j = 0; j < n; j++) ipiv(j) = 0;
    for (i=0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++)
            if (ipiv(j) != 1)
                for (k = 0; k < n; k++) {
                    if (ipiv(k) == 0) {
                        if (abs(a(j,k)) >= big) {
                            big = abs(a(j,k));
                            irow = j;
                            icol = k;
                        }
                    }
                }
        ++(ipiv(icol));
        if (irow != icol) {
            for (l = 0; l < n; l++) SWAP(a(irow, l), a(icol, l));
            SWAP(b(irow), b(icol));
        }
        indxr(i) = irow;
        indxc(i) = icol;
        if (a(icol, icol) == 0.0) throw("gaussj: Singular Matrix");
        pivinv = 1.0/a(icol, icol);
        a(icol, icol) = 1.0;
        for (l = 0; l < n; l++) a(icol, l) *= pivinv;
        b(icol) *= pivinv;
        for (ll = 0; ll < n; ll++)
            if (ll != icol) {
                dum = a(ll, icol);
                a(ll, icol) = 0.0;
                for (l = 0; l < n; l++) a(ll, l) -= a(icol, l)*dum;
                b(ll) -= b(icol)*dum;
            }
    }
    for (l = n - 1; l >= 0; l--) {
        if (indxr(l) != indxc(l))
            for (k = 0; k < n; k++)
                SWAP(a(k, indxr(l)), a(k, indxc(l)));
    }
    return;
}
