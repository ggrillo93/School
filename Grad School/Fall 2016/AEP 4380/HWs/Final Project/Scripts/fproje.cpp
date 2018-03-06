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
#include "nr3.h"
#define ARRAYT_BOUNDS_CHECK
using namespace std;

int main() {
	int N = 2, M = 200, i, j, k, n1 = 1, n2 = 1, nite = 0;
	double K = 2.478e14, totmass = 1.9891e30, G = 66732e-8, err = 0;
	arrayt<double> Evec(N*M), negEvec(N*M), mlist(M), y(M,N), deltay(M*N), relaxmat(M*N, M*N), sjnmatint(N,2*N), sjnmat0(n1, N), sjnmats(n2, N);
	ofstream fp, sp;
	void Sjn(int, arrayt<double>&, arrayt<double>&, arrayt<double>&, double);
	double Ecompint(int, int, arrayt<double>&, arrayt<double>&, double);
	double Ecomp0(arrayt<double>&, arrayt<double>&, double);
	double EcompM(arrayt<double>&, arrayt<double>&, double);
	void gaussj(arrayt<double>&, arrayt<double>&);
	for (k = 0; k < M; k++) {
		mlist(k) = k*totmass/M;
		y(k,0) = k*1e10/M;
	}
	for (k = 0; k < M; k++) { y(M-1-k,1) = k*1e20/M; }
	for (k = 1; k < M; k++) { // initialize interior points of Evec
		Evec(2*k-1) = Ecompint(0, k, y, mlist, K);
		Evec(2*k) = Ecompint(1, k, y, mlist, K);
	}
	Evec(0) = Ecomp0(y, mlist, K);
	Evec(N*M-1) = EcompM(y, mlist, K);
	// for (k = 0; k < M*N; k++) {cout << Evec(k) << endl; }
	for (k = 0; k < M; k++) {
		err += abs(Evec(2*k-1)/1e10);
		err += abs(Evec(2*k)/1e14);
	}
	err /= Evec.n();
	cout << err << endl;
	while (err > 1 && nite < 1000) {
		Sjn(0, y, mlist, sjnmat0, K);
		Sjn(M, y, mlist, sjnmats, K);
		relaxmat(0,0) = sjnmat0(0,0);
		relaxmat(0,1) = sjnmat0(0,1);
		relaxmat(M*N-1,M*N-2) = sjnmats(0,0);
		relaxmat(M*N-1,M*N-1) = sjnmats(0,1);
		for (k = 1; k < M; k++) {
			Sjn(k, y, mlist, sjnmatint, K);
			relaxmat(2*k-1, 2*k-2) = sjnmatint(0,0);
			relaxmat(2*k-1, 2*k-1) = sjnmatint(0,1);
			relaxmat(2*k-1, 2*k) = sjnmatint(0,2);
			relaxmat(2*k-1, 2*k+1) = sjnmatint(0,3);
			relaxmat(2*k, 2*k-2) = sjnmatint(1,0);
			relaxmat(2*k, 2*k-1) = sjnmatint(1,1);
			relaxmat(2*k, 2*k) = sjnmatint(1,2);
			relaxmat(2*k, 2*k+1) = sjnmatint(1,3);
		}
		for (k = 0; k < N*M; k++) { negEvec(k) = -Evec(k); }
		gaussj(relaxmat, negEvec);
		for (k = 0; k < M*N; k++) { deltay(k) = negEvec(k); }
		for (k = 0; k < M; k++) {
			y(k,0) = y(k,0) + deltay(2*k-1);
			y(k,1) = y(k,1) + deltay(2*k);
		}
		// for (k = 0; k < M; k++) {cout << y(k,0) << setw(15) << y(k,1) << setw(15) << mlist(k) << endl;}
		for (k = 1; k < M; k++) { // initialize interior points of Evec
			Evec(2*k-1) = Ecompint(0, k, y, mlist, K);
			Evec(2*k) = Ecompint(1, k, y, mlist, K);
		}
		Evec(0) = Ecomp0(y, mlist, K);
		Evec(N*M-1) = EcompM(y, mlist, K);
		// for (k = 0; k < M*N; k++) {cout << Evec(k) << endl;}
		err = 0;
		for (k = 0; k < M; k++) {
			err += abs(Evec(2*k-1)/1e10);
			err += abs(Evec(2*k)/1e14);
		}
		err /= Evec.n();
		cout << err << endl;
		nite++;
	}
	fp.open("solution.dat");
	for (k = 0; k < M; k++) {
		fp << y(k,0) << setw(15) << y(k,1) << setw(15) << mlist(k) << endl;
	}
	fp.close();
	return (EXIT_SUCCESS);
}

double Ecomp0(arrayt<double>& y, arrayt<double>& mlist, double K) {
	return 0.0;
}

double EcompM(arrayt<double>& y, arrayt<double>& mlist, double K) {
	return 0.0;
}

double Ecompint(int j, int k, arrayt<double>& y, arrayt<double>& mlist, double K) {
	double pi = 4.0 * atan(1.0), G = 6.6732e-8, pol = 3.0/5.0, ans;
	if (j == 0) {
		ans = -pow(2*K,pol)*(mlist(k)-mlist(k-1))/(pi*pow(y(k-1,1) + y(k,1), pol)*pow(y(k-1,0)+y(k,0),2.0))-y(k-1,0)+y(k,0);
		return ans;
	}
	else if (j == 1) {
		return 2*G*(mlist(k)-mlist(k-1))*(mlist(k-1)+mlist(k))/(pi*pow(y(k-1,0)+y(k,0),4.0)-y(k-1,1)+y(k,1));
	}
}

void Sjn(int k, arrayt<double>& y, arrayt<double>& mlist, arrayt<double>& sjnmat, double K) {
	int N = y.n2(), M = mlist.n();
	double G = 6.6732e-8, pi = 4.0 * atan(1.0), S00, S01, S10;
	if (k == 0) {
		sjnmat(0,0) = 1.0;
		sjnmat(0,1) = 0.0;
	}
	else if (k == M) {
		sjnmat(0,0) = 0.0;
		sjnmat(0,1) = 1.0;
	}
	else {
		S00 = pow(2.0, 8.0/5.0)*pow(K,3.0/5.0)*(mlist(k)-mlist(k-1))/ (pi*pow(y(k-1,1)+y(k,1),3.0/5.0)*pow(y(k-1,0)+y(k,0),3.0))-1;
		S01 = 3.0*pow(2.0,3.0/5.0)*pow(K,3.0/5.0)*(mlist(k)-mlist(k-1))/ (5*pi*pow(y(k-1,1)+y(k,1),8.0/5.0)*pow(y(k-1,0)+y(k,0),2.0));
		S10 = -8.0*G*(mlist(k)-mlist(k-1))*(mlist(k-1)*mlist(k))/(pi*pow(y(k-1,0)+y(k,0),5.0));
		sjnmat(0,0) = S00;
		sjnmat(0,1) = S01;
		sjnmat(0,2) = S00+2;
		sjnmat(0,3) = S01;
		sjnmat(1,0) = S10;
		sjnmat(1,1) = -1.0;
		sjnmat(1,2) = S10;
		sjnmat(1,3) = 1.0;
	}
	return;
}

void gaussj(arrayt<double> &a, arrayt<double> &b)
{
	int i,icol,irow,j,k,l,ll,n=a.n1(),m=b.n2();
	double big,dum,pivinv;
	arrayt<int> indxc(n),indxr(n),ipiv(n);
	for (j=0;j<n;j++) ipiv(j)=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv(j) != 1)
				for (k=0;k<n;k++) {
					if (ipiv(k) == 0) {
						if (abs(a(j,k)) >= big) {
							big=abs(a(j,k));
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv(icol));
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a(irow,l),a(icol,l));
			for (l=0;l<m;l++) SWAP(b(irow,l),b(icol,l));
		}
		indxr(i)=irow;
		indxc(i)=icol;
		if (a(icol,icol) == 0.0) throw("gaussj: Singular Matrix");
		pivinv=1.0/a(icol,icol);
		a(icol,icol)=1.0;
		for (l=0;l<n;l++) a(icol,l) *= pivinv;
		for (l=0;l<m;l++) b(icol,l) *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a(ll,icol);
				a(ll,icol)=0.0;
				for (l=0;l<n;l++) a(ll,l) -= a(icol,l)*dum;
				for (l=0;l<m;l++) b(ll,l) -= b(icol,l)*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr(l) != indxc(l))
			for (k=0;k<n;k++)
				SWAP(a(k,indxr(l)),a(k,indxc(l)));
	}
	return;
}
//
// void mnewt(const int ntrial, arrayt<double> &x, const double tolx, const double tolf, arrayt<double> x, arrayt<double> fvec, arrayt<double> fjac) {
// 	int i,n=x.n2();
// 	arrayt<double> p(n);
// 	for (Int k=0;k<ntrial;k++) {
// 		double errf=0.0;
// 		for (i=0;i<n;i++) errf += abs(fvec(i));
// 		if (errf <= tolf) return;
// 		for (i=0;i<n;i++) p(i) = -fvec(i);
// 		LUdcmp alu(fjac);
// 		alu.solve(p,p);
// 		double errx=0.0;
// 		for (i=0;i<n;i++) {
// 			errx += abs(p(i));
// 			x(i) += p(i);
// 		}
// 		if (errx <= tolx) return;
// 	}
// 	return;
// }
