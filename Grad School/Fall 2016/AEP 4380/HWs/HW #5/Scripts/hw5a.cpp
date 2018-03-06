/* AEP 4380 HW#5a
  Autostep size RK and Planetary Motion
  Run on core i7 with g++ 5.4.0 (Ubuntu)
  Gianfranco Grillo 10/14/2016
*/

#include <cstdlib> // plain C
#include <cmath>

#include <iostream> // stream IO
#include <fstream> // stream file IO
#include <iomanip> // to format the output

using namespace std;

int main() {
	int istep, n = 16, i;
	double initvar[n] = {149.598e9, 228.0e9, 778.29e9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 2.9786e4, 2.4127e4, 1.30588e4, -3.0e1};
	double finalvar[16], t0 = 0.0, hold, h, nsecinday, totalsec, nsecinyear;
	void autosteprk4(double[], double[], int, double&, double, void(double[], double, double[], double));
	void nbody(double[], double, double[], double);
	ofstream fp;
    fp.open( "solarsystem1.dat" ); // open new file for output
	if( fp.fail() ) { // in case of error
        cout << "cannot open file" << endl;
        return( EXIT_SUCCESS );
	}
	nsecinyear = 365.0*24.0*60.0*60.0;
	nsecinday = 3600.0*24.0;
	h = nsecinday; // initial h
	totalsec = 15.0*nsecinyear;
	fp << setw(15) << "t (years)" << setw(15) << "xearth" << setw(15) << "yearth" << setw(15) << "xmars" << setw(15)
	<< "ymars" << setw(15) << "xjupiter" << setw(15) << "yjupiter" << setw(15) << "xsun" << setw(15) << "ysun" << endl;
	while (t0 < totalsec) {
		hold = h;
		autosteprk4(initvar, finalvar, n, h, t0, nbody);
		if (h > hold || h == hold) { // means error is low enough, result is correct
			for (i = 0; i < n; i++) { initvar[i] = finalvar[i]; }
			fp <<  setw(15) <<  t0/nsecinyear << setw(15) << finalvar[0] << setw(15) << finalvar[4] << setw(15)
			<< finalvar[1] << setw(15) << finalvar[5] << setw(15) << finalvar[2] << setw(15) << finalvar[6]
			<< setw(15) << finalvar[3] << setw(15) << finalvar[7] << endl;
			t0 = t0 + hold;
		}
  }
    fp.close();
    return( EXIT_SUCCESS );
}

void autosteprk4(double initvar[], double highfinalvar[], int n, double& h,
				 double t0, void functions(double[], double, double[], double)) {
	int i;
	double scale[n], error[n], maxerror = 5e-8, lowfinalvar[n], deltamax;
	double maximum(double[], int);
	double *k1, *k2, *k3, *k4, *k5, *k6, *k7, *scalek, *temp;
	double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0/9.0, c6 = 1.0, c7 = 1.0;
	double a21 = 0.2, a31 = 3.0/40.0, a32 = 9.0/40.0, a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0;
	double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0;
	double a61 = 9017.0/3168.0, a62 = -355.0/33.0, a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0;
	double b1 = 35.0/384.0, b3 = 500.0/1113.0, b4 = 125.0/192.0, b5 = -2187.0/6784.0, b6 = 11.0/84.0;
	double d1 = 5179.0/57600.0, d3 = 7571.0/16695.0, d4 = 393.0/640.0, d5 = -92097.0/339200.0, d6 = 187.0/2100.0, d7 = 1.0/40.0;
	k1 = new double[10*n];
	if (NULL == k1) {
		cout << "can't allocate arrays in rk4" << endl;
	return;
	}
	k2 = k1 + n;
	k3 = k2 + n;
	k4 = k3 + n;
	k5 = k4 + n;
	k6 = k5 + n;
	k7 = k6 + n;
	temp = k7 + n;
	scalek = temp + n; // allocate dynamic memory to variables
	functions(initvar, t0, k1, h); // generate k1
	for (i = 0; i < n; i++) {
		temp[i] = initvar[i] + a21*k1[i];
	}
	functions(temp, t0 + c2*h, k2, h);
	for (i = 0; i < n; i++) {
		temp[i] = initvar[i] + a31*k1[i] + a32*k2[i];
	}
	functions(temp, t0 + c3*h, k3, h);
	for (i = 0; i < n; i++) {
		temp[i] = initvar[i] + a41*k1[i] + a42*k2[i] + a43*k3[i];
	}
	functions(temp, t0 + c4*h, k4, h);
	for (i = 0; i < n; i++) {
		temp[i] = initvar[i] + a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i];
	}
	functions(temp, t0 + c5*h, k5, h);
	for (i = 0; i < n; i++) {
		temp[i] = initvar[i] + a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i];
	}
	functions(temp, t0 + c6*h, k6, h);
	for (i = 0; i < n; i++) {
		highfinalvar[i] = initvar[i] + b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i];
	}
	functions(highfinalvar, t0 + h, k7, h);
	for (i = 0; i < n; i++) {
		lowfinalvar[i] = initvar[i] + d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k7[i];
	}
	functions(initvar, t0 + h, scalek, h); // generate scale
	for (i = 0; i < n; i++) {
		scale[i] = abs(initvar[i]) + abs(scalek[i]) + 0.01;
	}
	for (i = 0; i < n; i++) {
		error[i] = abs((highfinalvar[i] - lowfinalvar[i])/scale[i]);
	}
	deltamax = maximum(error, n);
	if (5*deltamax < maxerror) {
		h = h*pow(maxerror/deltamax, 0.2); // if error is too small, increase h
		delete k1;
		return;
	}
	else if (deltamax > maxerror && h > 1e-12) { // if error is too big, decrease h
		h = h/3.0;
		delete k1;
		return;
	}
	else {
		delete k1;
		return;
	}
}

void nbody(double var[], double t, double k[16], double h) {
	int nobj = 4, i, n, vxloc, vyloc;
	double masses[nobj] = {5.9742e24, 0.64191e24, 1898.8e24, 1.9891e30};
	double G = 6.6726e-11, xsingleacc, xdistance, ysingleacc, ydistance, totaldist, totaldistc;
	for (i = 0; i < nobj; i++) {
		vxloc = i+2*nobj;
		vyloc = i+3*nobj;
		k[i] = h*var[i+2*nobj];
		k[i+nobj] = h*var[i+3*nobj];
		k[vxloc] = 0;
		k[vyloc] = 0;
		for(n = 0; n < nobj; n++) {
			if (n != i) {
				xdistance = var[n] - var[i];
				ydistance = var[n+nobj] - var[i+nobj];
				totaldist = sqrt(xdistance*xdistance+ydistance*ydistance);
				totaldistc = totaldist*totaldist*totaldist;
				xsingleacc = h*G*masses[n]*xdistance/totaldistc;
				k[vxloc] = k[vxloc] + xsingleacc;
				ysingleacc = h*G*masses[n]*ydistance/totaldistc;
				k[vyloc] = k[vyloc] + ysingleacc;
			}
		}
	}
	return;
}

double maximum(double array[], int n) {
	double max;
	int i;
	max = array[0];
	for (i = 0; i < n; i++) {
		if (array[i] > max) {
			max = array[i];
		}
	}
	return max;
}
