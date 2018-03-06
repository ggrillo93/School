/* AEP 4380 HW#6f
  Boundary value problems and relaxation
  Run on core i7 with g++ 5.4.0 (Ubuntu)
  Gianfranco Grillo 10/25/2016
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <iostream> // stream IO
#include <fstream> // stream file IO
#include <iomanip> // to format the output
#define ARRAYT_BOUNDS_CHECK
#include "arrayt.hpp"

using namespace std;

int main() {
    float rmax = 15.0, zmax = 30.0, z, r, zlength, rlength, maxdiff = 2000.0, w = 1.84, h = 0.06, tol = 1.5, Er, Ez, Esqr, eps = 8.8543e-3, totalEsqr, C1, C2;
    int i, j, nite, t;
    zlength = zmax/h;
    rlength = rmax/h;
    arrayt<float> oldgrid(zlength, rlength);
    arrayt<char> elecloc(zlength, rlength);
    void relax(arrayt<float>&, arrayt<char>&, float, float&);
    ofstream fp, rp, zp;
    rp.open("cplotr.dat");
    zp.open("cplotz.dat");
    fp.open( "cplotmain.dat" );
    for (i = 0; i < zlength; i++) for (j = 0; j < rlength; j++) { // set initial values of potential and flag electrodes
        z = i*h;
        r = j*h;
        if (z == 13.0 && r >= 2.0 && r <= 10.0) {
            elecloc(i,j) = 1;
            oldgrid(i,j) = 0.0;
        }    
        else if (z == 15.0 && r >= 3.0 && r <= 10.0) {
            elecloc(i,j) = 1;
            oldgrid(i,j) = 2000.0;
        }
        else if (z >= 16.5 && z <= 18 && r >= 3.0 && r <= 10.0) {
            elecloc(i,j) = 1;
            oldgrid(i,j) = 0.0;
        } else {
            elecloc(i,j) = 0;
            oldgrid(i,j) = 0.0;
        }
    }
    nite = 0;
    maxdiff = 2000.0;
    while ((maxdiff > tol) && (nite < 10000)) {
        if (maxdiff > 2000.1) {
            nite = 10000;
            break;
        }
        relax(oldgrid, elecloc, w, maxdiff);
        nite++;
    }
    totalEsqr = 0.0;
    for (i = (13.0+h)/h; i < 15.0/h; i++) for (j=1; j < rlength-1; j++) {
        z = i*h;
        r = j*h;
        Ez = (oldgrid(i+1,j)-oldgrid(i-1,j))/(2*h);
        Er = (oldgrid(i,j+1)-oldgrid(i,j-1))/(2*h);
        Esqr = Ez*Ez+Er*Er;
        totalEsqr = totalEsqr + Esqr;
    }
    C1 = (eps*totalEsqr*2*3.1415)/(2000*2000);
    cout << C1 << endl;
    totalEsqr = 0.0;
    for (i = (15.0+h)/h; i < 16.5/h; i++) for (j=1; j < rlength-1; j++) {
        z = i*h;
        r = j*h;
        Ez = (oldgrid(i+1,j)-oldgrid(i-1,j))/(2*h);
        Er = (oldgrid(i,j+1)-oldgrid(i,j-1))/(2*h);
        Esqr = Ez*Ez+Er*Er;
        totalEsqr = totalEsqr + Esqr;
    }
    C2 = (eps*totalEsqr*2*3.1415)/(2000*2000);
    cout << C2 << endl;
    return( EXIT_SUCCESS );
}

void relax(arrayt<float>& oldgrid, arrayt<char> &flags, float w, float& maxdiff) {
    int zlength = oldgrid.n1(), rlength = oldgrid.n2(), i, j, n;
    float max, diff;
    arrayt<float> fdgrid(zlength, rlength);
    max = 0.0;
    for (i = 1; i < zlength-1; i++) {
        for (j = 0; j < rlength-1; j++) {
            if (flags(i,j) == 1) {
                fdgrid(i,j) = oldgrid(i,j);
            }
            else if (j == 0) {
                fdgrid(i,j) = (4*oldgrid(i,1) + oldgrid(i+1,0) + oldgrid(i-1,0))/6.0;
            } else {
                fdgrid(i,j) = (0.25*(oldgrid(i,j+1) + oldgrid(i,j-1) + oldgrid(i+1,j) +
                                             oldgrid(i-1,j)) + (oldgrid(i,j+1) - oldgrid(i,j-1))/(8*j));
            }
            oldgrid(i,j) = oldgrid(i,j) + w*(fdgrid(i,j)-oldgrid(i,j));
            diff = abs(fdgrid(i,j)-oldgrid(i,j));
            oldgrid(i,j) = oldgrid(i,j) + w*(fdgrid(i,j)-oldgrid(i,j));
            if (diff > max) {	max = diff; }
        }
    }
    maxdiff = max;
    return;
}
