/* AEP 4380 HW#3a
  Bessel Function Plotting
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
    double x, J0, J1, J2, Y0, Y1, Y2;
    ofstream fp;
    fp.open( "jfunctions.dat" ); // open new file for output
	if( fp.fail() ) { // in case of error
        cout << "cannot open file" << endl;
        return( EXIT_SUCCESS );
    }
    for (x=0.002; x < 20.0; x = x + 0.002) { // loop between 0 and 20, get 10000 points
        J0 = myBessel.j0(x);
        J1 = myBessel.j1(x);
        J2 = myBessel.jn(2,x);
        fp << setw(15) << x << setw(15) << J0 << setw(15) << J1 << setw(15) << J2 << endl; // output columns to file
    }
    fp.close();
    
    fp.open( "yfunctions.dat" );
    if( fp.fail() ) { // in case of error
        cout << "cannot open file" << endl;
        return( EXIT_SUCCESS );
    }
    for (x=0.76925; x < 20; x = x + 0.001925) { // loop between 0.75 and 20, get 10000 points
        Y0 = myBessel.y0(x);
        Y1 = myBessel.y1(x);
        Y2 = myBessel.yn(2,x);
        fp << setw(15) << x << setw(15) << Y0 << setw(15) << Y1 << setw(15) << Y2 << endl;
    }
    fp.close();
    return( EXIT_SUCCESS );
 }