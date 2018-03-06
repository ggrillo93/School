/* 	AEP 438 Example #1
	Test numerical derivatives
*/

#include <cstdlib> // plain C
#include <cmath>

#include <iostream> // stream IO
#include <fstream> // stream file IO
#include <iomanip> // to format the output

using namespace std;

int main()
{
	int i, n=200;
	double h=0.5,xmin=-7.0, xmax=+7.0, x, dx, f1, fpbd, fpfd, fpcd;
	double feval( double );
	ofstream fp;		// output file using streams
	fp.open( "backwarddiff.dat"); // open new file for output
	if( fp.fail() ) { // or fp.bad()
		cout << "cannot open file" << endl;
		return( EXIT_SUCCESS );
	}
	
	dx = ( xmax - xmin ) / (n-1);
	for( i=0; i<n; i++ ){
		x = xmin + i * dx;
		f1 = feval(x);
		fpbd = ( f1 - feval(x-h) )/h;
		
		// data file for python, Matlab; need to separate into col.
		fp << setw(15) << x << setw(15) << f1 << setw(15) << fpbd << endl;
		
		// print to screen
		fp << setw(15) << x << setw(15) << f1 << setw(15) << fpbd << endl;
	}
	
	fp.close();
	
	fp.open("forwarddiff.dat");
		if( fp.fail() ) { // or fp.bad()
		cout << "cannot open file" << endl;
		return( EXIT_SUCCESS );
	}
	
	dx = ( xmax - xmin ) / (n-1);
	for( i=0; i<n; i++ ){
		x = xmin + i * dx;
		f1 = feval(x);
		fpfd = ( feval(x+h) - f1 )/h;
		
		// data file for python, Matlab; need to separate into col.
		fp << setw(15) << x << setw(15) << f1 << setw(15) << fpfd << endl;
		
		// print to screen
		fp << setw(15) << x << setw(15) << f1 << setw(15) << fpfd << endl;
	}
	
	fp.close();
	
	fp.open("centraldiff.dat");
		if( fp.fail() ) { // or fp.bad()
		cout << "cannot open file" << endl;
		return( EXIT_SUCCESS );
	}
	
	dx = ( xmax - xmin ) / (n-1);
	for( i=0; i<n; i++ ){
		x = xmin + i * dx;
		f1 = feval(x);
		fpcd = ( feval(x+h) - feval(x-h))/(2*h);
		
		// data file for python, Matlab; need to separate into col.
		fp << setw(15) << x << setw(15) << f1 << setw(15) << fpcd << endl;
		
		// print to screen
		fp << setw(15) << x << setw(15) << f1 << setw(15) << fpcd << endl;
	}
	
	fp.close();
	
	return( EXIT_SUCCESS );

} // end main()

double feval( double x )
{ return( sin(x) * exp( -0.04*x*x ) ); }