/* AEP 4380 HW#6a
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
	float h = 0.1, rmax = 15.0, zmax = 30.0, z, r, zlength, rlength, tol = 2.0, maxdiff;
	int i, j, nite;
	double w;
	void relax(arrayt<float>&, arrayt<char>&, float, float&);
	zlength = zmax/h;
	rlength = rmax/h;
	arrayt<float> oldgridc(zlength, rlength);
	arrayt<float> oldgrid(zlength, rlength);
	arrayt<char> elecloc(zlength, rlength);
	ofstream fp;
    fp.open( "wvsiterations2.dat" );
	for (i = 0; i < zlength; i++) {
		for (j = 0; j < rlength; j++) { // set initial values of potential and flag electrodes
			z = i*h;
			r = j*h;
			if (z == 13.0 && r >= 2.0 && r <= 10.0) {
				elecloc(i,j) = 1;
				oldgridc(i,j) = 0.0;
			}	
			else if (z == 15.0 && r >= 3.0 && r <= 10.0) {
				elecloc(i,j) = 1;
				oldgridc(i,j) = 2000.0;
			}
			else if (z >= 16.5 && z <= 18 && r >= 3.0 && r <= 10.0) {
				elecloc(i,j) = 1;
				oldgridc(i,j) = 0.0;
			}
			else {
				elecloc(i,j) = 0;
				oldgridc(i,j) = 0.0;
			}
		}
	}
	for (w = 1.0; w < 2.0; w = w + 0.01) {
		nite = 0;
		maxdiff = 2000.0;
		oldgrid = oldgridc;
		while ((maxdiff > tol) && (nite < 10000)) {
			if (maxdiff > 2000.0) {
				nite = 10000;
				break;
			}
			relax(oldgrid, elecloc, w, maxdiff);
			nite++;
		}
		fp << setw(15) << w << setw(15) << nite << endl;
	}
    fp.close();
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
				fdgrid(i,j) = ((4*oldgrid(i,1) + oldgrid(i+1,0) + oldgrid(i-1,0))/6.0);
			}
			else {
				fdgrid(i,j) = (0.25*(oldgrid(i,j+1) + oldgrid(i,j-1) + oldgrid(i+1,j) + oldgrid(i-1,j)) + (oldgrid(i,j+1) - oldgrid(i,j-1))/(8*j));
			}
			diff = abs(fdgrid(i,j)-oldgrid(i,j));
			oldgrid(i,j) = oldgrid(i,j) + w*(fdgrid(i,j)-oldgrid(i,j));
			if (diff > max) {
				max = diff;
			}
		}
	}
	maxdiff = max;
	return;
}
