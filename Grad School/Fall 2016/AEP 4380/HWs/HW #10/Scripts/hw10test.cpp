/* AEP 4380 HW#10a
  Fast Fourier Transform and Spectral Methods
  Run on core i7 with g++ 5.4.0 (Ubuntu)
  Gianfranco Grillo 11/26/2016
*/

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <fftw3.h>
#include "arrayt.hpp"

using namespace std;

int main() {
    fftw_complex *psi;
    fftw_plan planTf, planTi;
    int Nx = 512, Ny = 512, Lx = 500, Ly = 500, i, j;
    arrayt<double> psilist(Nx);
    double x, y, arg1, arg2, arg3, x1n, y1n, x2n, y2n, x3n;
    psi = (fftw_complex*) fftw_malloc(Nx*sizeof(fftw_complex));
    ofstream fp;

    planTf = fftw_plan_dft_1d(Nx, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE); // forward
    planTi = fftw_plan_dft_1d(Nx, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE); // inverse

    for (i = 0; i < Nx; i++) for (j = 0; j < Ny; j++) {
        x = i*Lx/Nx;
        psi[i][0] = sin(200*x); // initial real part
        psi[i][1] = 0; // initial imaginary part
    }
    fftw_execute_dft(planTf, psi, psi);
    fp.open("2ndpsitest.dat");
    for (i = 0; i < Nx; i++) {
        fp << psi[i][1] << endl;
    }
    fp.close();
    return(EXIT_SUCCESS);
}
