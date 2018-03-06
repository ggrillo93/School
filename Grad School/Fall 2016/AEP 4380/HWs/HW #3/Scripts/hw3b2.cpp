/* AEP 4380 HW#3b2
  Root finding of Bessel Functions
  Run on core i7 with g++ 4.8.1
  Gianfranco Grillo 09/19/2016
*/

#include <cstdlib> // plain C
#include <cmath>

#include <iostream> // stream IO
#include <fstream> // stream file IO
#include <iomanip> // to format the output
#include "nr3.h"
#include "bessel.h"

using namespace std;
Bessjy myBessel; // make instance of the Bessel function

int main()
{
    // calculate roots using bisection
    double x11=0.1, x12=1.0, x21=2.5, x22=3.5, x31=4.0; // can't use x11=0 because 0 has no sign...
    double x32=5.0, x41=6.0, x42=7.0, x51=7.0, x52=8.0;
    double b1, b2, b3, b4, b5, tol;
    int it1 = 0, it2 = 0; // set iteration count
    double g(double);
    double bisect(double(*)(double), double, double, double, int&);
    ofstream fp;

    tol = 0.000001; // set tolerance value
    
    b1 = bisect(g, x11, x12, tol, it1);
    b2 = bisect(g, x21, x22, tol, it1);
    b3 = bisect(g, x31, x32, tol, it1);
    b4 = bisect(g, x41, x42, tol, it1);
    b5 = bisect(g, x51, x52, tol, it1);
    
    // calculate roots using false position
    double r1, r2, r3, r4, r5;
    double rf(double(*)(double), double, double, double, int&);
    r1 = rf(g, x11, x12, tol, it2);
    r2 = rf(g, x21, x22, tol, it2);
    r3 = rf(g, x31, x32, tol, it2);
    r4 = rf(g, x41, x42, tol, it2);
    r5 = rf(g, x51, x52, tol, it2);
    
    // output to file
    fp.open( "roots.dat");
	if( fp.fail() ) {
        cout << "cannot open file" << endl;
		return( EXIT_SUCCESS );
	}
    fp << setw(15) << "Bisection method" << setw(15) << b1 << setw(15) << b2 << setw(15) << b3 << setw(15) << b4 << setw(15) << b5 << setw(15) << it1 << endl;
    fp << setw(15) << "False position method" << setw(15) << r1 << setw(15) << r2 << setw(15) << r3 << setw(15) << r4 << setw(15) << r5 << setw(15) << it2<< endl;
    fp.close();
    return( EXIT_SUCCESS );
 } // end main
 
 double g(double x) { // calculates J0(x)Y0(x)-J2(x)Y2(x)
    double J0, Y0, J2, Y2;
    J0 = myBessel.j0(x);
    J2 = myBessel.jn(2,x);
    Y0 = myBessel.y0(x);
    Y2 = myBessel.yn(2,x);
    return J0*Y0-J2*Y2;
 }
 
 double bisect(double(*g)(double), double x1, double x2, double tol, int& nit) {
    double x3, g1, g2, g3;
    double bisect(double(*)(double), double, double, double, int&);
    nit++; // increase iteration count
    g1 = g(x1);
    g2 = g(x2);
    if((g1 > 0 && g2 > 0) || (g1 < 0 && g2 < 0)) { // check whether g1, g2 have opposite signs
        cout << "No root in interval x1 < x < x2";
        return( EXIT_SUCCESS );
    }
    x3 = 0.5*(x1+x2); // calculate x3 by averaging x1 and x2
    g3 = g(x3);
    cout << setw(15) << x3 << setw(15) << g3 << endl; // for verification purposes
    if(abs(g3) > tol) { // check whether to continue looping
        if((g3 > 0 && g1 > 0) || (g3 < 0 && g1 < 0)) {
            x1 = x3;
            bisect(g,x1,x2,tol, nit);
        }
        else if((g3 > 0 && g2 > 0) || (g3 < 0 && g2 < 0)) {
            x2 = x3;
            bisect(g,x1,x2,tol, nit);
        }
    }
    else {
        return x3; // once looping is done, return root
    }
 }
 
double rf(double(*g)(double), double x1, double x2, double tol, int& nit) {
    double x3, g1, g2, g3;
    double rf(double(*)(double), double, double, double, int&);
    nit++; // increase iteration count
    g1 = g(x1);
    g2 = g(x2);
    if((g1 > 0 && g2 > 0) || (g1 < 0 && g2 < 0)) { // check whether g1, g2 have opposite signs
        cout << "No root in interval x1 < x < x2";
        return( EXIT_SUCCESS );
    }
    x3 = x1 - g1*(x2-x1)/(g2-g1); // interpolate linearly
    g3 = g(x3);
    cout << setw(15) << x3 << setw(15) << g3 << endl; // for verification purposes
    if(abs(g3) > tol) { // check whether to continue looping
        if((g3 > 0 && g1 > 0) || (g3 < 0 && g1 < 0)) {
            x1 = x3;
            bisect(g,x1,x2,tol, nit);
        }
        else if((g3 > 0 && g2 > 0) || (g3 < 0 && g2 < 0)) {
            x2 = x3;
            bisect(g,x1,x2,tol, nit);
        }
    }
    else {
        return x3; // once looping is done, return root
    }
 }