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
    fftw_complex *psi, *coeff;
    fftw_plan planTf, planTi;
    int Nx = 512, Ny = 512, Lx = 500, Ly = 500, i, j;
    arrayt<double> psimat(Nx,Ny), k(Nx,Ny);
    arrayt<int> auxindex(Nx);
    double x, y, arg1, arg2, arg3, x1n, y1n, x2n, y2n, x3n, v = 343.0, kx, ky, pi = 4.0*atan(1.0), Nsqr = Nx*Ny;
    psi = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex));
    coeff = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex));
    ofstream fp;
    fp.open("t=0.dat");

    planTf = fftw_plan_dft_2d(Nx, Ny, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE); // forward
    planTi = fftw_plan_dft_2d(Nx, Ny, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE); // inverse

    for (i = 0; i < Nx; i++) for (j = 0; j < Ny; j++) {
        x = i*Lx/Nx;
        y = j*Ly/Ny;
        x1n = x-0.4*Lx;
        y1n = y-0.4*Ly;
        arg1 = -(x1n*x1n+y1n*y1n)/100.0;
        x2n = x-0.5*Lx;
        y2n = y-0.6*Ly;
        arg2 = -(x2n*x2n+y2n*y2n)/400.0;
        x3n = x-0.6*Lx;
        arg3 = -(x3n*x3n+y1n*y1n)/100.0;
        psi[i*Nx + j][0] = exp(arg1)+2*exp(arg2)+exp(arg3); // initial real part
        psimat(i,j) = psi[i*Nx + j][0];
        psi[i*Nx + j][1] = 0; // initial imaginary part
    }

    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) { fp << psimat(i, j) << "    "; }
        fp << endl;
    }
    fp.close();
    fftw_execute_dft(planTf, psi, psi);

    for (i = 0; i <= Nx/2-1; i++) { auxindex(i) = i; }
    j = 0;
    for (i = Nx/2; i <= Nx; i++) {
        auxindex(i) = -Nx/2 + j;
        j++;
    }
    for (i = 0; i < Nx; i++) for (j = 0; j < Ny; j++) {
        coeff[i*Nx + j][0] = psi[i*Nx + j][0];
        coeff[i*Nx + j][1] = psi[i*Nx + j][1];
        kx = 2*pi*auxindex(i)/Lx;
        ky = 2*pi*auxindex(j)/Ly;
        k(i,j) = sqrt(kx*kx+ky*ky);
        psi[i*Nx + j][0] *= cos(v*0.1*k(i,j));
        psi[i*Nx + j][1] *= cos(v*0.1*k(i,j));
    }

    fftw_execute_dft(planTi, psi, psi);
    fp.open("t=0.1.dat");
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            psimat(i,j) = psi[i*Nx + j][0]/Nsqr;
            fp << psimat(i, j) << "    ";
        }
        fp << endl;
    }
    fp.close();

    for (i = 0; i < Nx; i++) for (j = 0; j < Ny; j++) {
        psi[i*Nx + j][0] = coeff[i*Nx + j][0]*cos(v*0.2*k(i,j));
        psi[i*Nx + j][1] = coeff[i*Nx + j][1]*cos(v*0.2*k(i,j));
    }

    fftw_execute_dft(planTi, psi, psi);
    fp.open("t=0.2.dat");
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            psimat(i,j) = psi[i*Nx + j][0]/Nsqr;
            fp << psimat(i, j) << "    ";
        }
        fp << endl;
    }
    fp.close();

    for (i = 0; i < Nx; i++) for (j = 0; j < Ny; j++) {
        psi[i*Nx + j][0] = coeff[i*Nx + j][0]*cos(v*0.4*k(i,j));
        psi[i*Nx + j][1] = coeff[i*Nx + j][1]*cos(v*0.4*k(i,j));
    }

    fftw_execute_dft(planTi, psi, psi);
    fp.open("t=0.4.dat");
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            psimat(i,j) = psi[i*Nx + j][0]/Nsqr;
            fp << psimat(i, j) << "    ";
        }
        fp << endl;
    }
    fp.close();
    return(EXIT_SUCCESS);
}
