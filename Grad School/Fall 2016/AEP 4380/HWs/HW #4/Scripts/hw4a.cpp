/* AEP 4380 HW#4a
  Runge-Kutta and Chaos
  Run on core i7 with g++ 5.4.0 (Ubuntu)
  Gianfranco Grillo 10/03/2016
*/

#include <cstdlib> // plain C
#include <cmath>

#include <iostream> // stream IO
#include <fstream> // stream file IO
#include <iomanip> // to format the output

using namespace std;

int main() {
	int istep, nstep=30000, n = 2, i;
	double initvar[2] = {1.0,0.0}, finalvar[2], t, t0 = 0.0, h=0.001;
	void rk4(double[], double[], int, double, double, void(double[], double, double[], double));
	void harmonicf(double[], double, double[], double);
	ofstream fp;
    fp.open( "oscillator.dat" ); // open new file for output
	if( fp.fail() ) { // in case of error
        cout << "cannot open file" << endl;
        return( EXIT_SUCCESS );
	}
	for (istep = 0; istep < nstep; istep++) {
		t = t0 + h*istep;
		rk4(initvar, finalvar, n, h, t, harmonicf);
		for (i = 0; i < n; i++) { initvar[i] = finalvar[i]; }
        fp <<  setw(15) <<  t << setw(15) << finalvar[0] << setw(15) << finalvar[1] << setw(15) << endl;
    }
    fp.close();
    return( EXIT_SUCCESS );
}

void rk4(double initvar[], double finalvar[], int n, double h, double t0, void functions(double[], double, double[], double)) {
	int i;
	double *k0, *k1, *k2, *k3, *temp;
	k0 = new double[5*n];
	if (NULL == k0) {
		cout << "can't allocate arrays in rk4" << endl;
		return;
	}
	k1 = k0 + n;
	k2 = k1 + n;
	k3 = k2 + n;
	temp = k3 + n;
	functions(initvar, t0, k0, h);
	for (i = 0; i < n; i++) {
		temp[i] = initvar[i] + 0.5*k0[i];
	}
	functions(temp, t0 + 0.5*h, k1, h);
	for (i = 0; i < n; i++) {
		temp[i] = initvar[i] + 0.5*k1[i];
	}
	functions(temp, t0 + 0.5*h, k2, h);
	for (i = 0; i < n; i++) {
		temp[i] = initvar[i] + k2[i];
	}
	functions(temp, t0 + h, k3, h);
	for (i = 0; i < n; i++) {
		finalvar[i] = initvar[i] + (1.0/6.0)*(k0[i]+2*k1[i]+2*k2[i]+k3[i]);
	}
	delete k0;
	return;
}

void harmonicf(double var[], double t, double k[2], double h) {
	k[0] = h*var[1];
	k[1] = -h*var[0];
	return;
}