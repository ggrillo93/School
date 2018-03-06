/* AEP 4380 HW#3b1
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
    double x, J0, J2, Y0, Y2;
    ofstream fp;
    fp.open( "rootplot.dat" ); // open new file for output
	if( fp.fail() ) { // in case of error
        cout << "cannot open file" << endl;
        return( EXIT_SUCCESS );
    }
    for (x=0.0; x < 20.0; x = x + 0.002) { // loop between 0 and 20, get 10000 points
        J0 = myBessel.j0(x);
        J2 = myBessel.jn(2,x);
        Y0 = myBessel.y0(x);
        Y2 = myBessel.yn(2,x);
        fp << setw(15) << x << setw(15) << J0*Y0 << setw(15) << J2*Y2 << endl; // output columns to file
    }
    fp.close();
    return( EXIT_SUCCESS );
 }