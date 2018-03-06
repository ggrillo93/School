/* AEP 4380 HW#2a
  Numerical Integration of Fresnel Integrals
  Run on core i7 with g++ 4.8.1
  Gianfranco Grillo 09/12/2016
*/

 #include <cstdlib>
 #include <cmath>
 
 #include <iostream>
 #include <fstream>
 #include <iomanip>

 using namespace std;
 
 int main() {
    double u0, I;
    int n;
    cout << "u0 = "; // ask for user input
    cin >> u0;
    double iint(double,int);
    ofstream fp;
	 fp.open( "fresnelintegral.dat"); // open new file for output
	 if( fp.fail() ) { // in case of error
		cout << "cannot open file" << endl;
		return( EXIT_SUCCESS );
	}
    for (n = 4; n < 8193; n = n*2) { // loop until n = 8192
      I = iint(u0,n);
      fp << setw(15) << n << setw(15) << I << endl; // output columns to file
    }
    
    fp.close();
	
	return( EXIT_SUCCESS );
 }
 
 double iint(double u0, int n) { // function that calculates I/I0
   double cint(double,int), C;
   double sint(double,int), S;
   C = cint(u0,n); // call functions that calculate C(u0), S(u0)
   S = sint(u0,n);
   return 0.5*(pow((C+0.5),2.0)+pow((S+0.5),2.0));
 }
 
 double cint(double u0, int n) { // calculates C(u0) with n points using trapezoidal rule
    double sum = 0.0, i;
    double sqrun(int,int,double), sq1, sq2;
    double halfpi = 2.0 * atan(1.0), ci;
    for (i = 1; i < n; i++) {
      sq1 = sqrun(i,n,u0);
      sq2 = sqrun(i+1,n,u0);
      ci = 0.5*(cos(halfpi*sq1)+cos(halfpi*sq2))*u0/(n-1);
      sum=sum+ci;
    }
    return sum;
 }
 
 double sint(double u0, int n) { // calculates S(u0) with n points using trapezoidal rule
    double sum = 0.0, i;
    double sqrun(int,int,double), sq1, sq2;
    double halfpi = 2.0 * atan(1.0), si;
    for (i = 1; i < n; i++) {
      sq1 = sqrun(i,n,u0);
      sq2 = sqrun(i+1,n,u0);
      si = 0.5*(sin(halfpi*sq1)+sin(halfpi*sq2))*u0/(n-1);
      sum=sum+si;
    }
    return sum;
 }
 
 double sqrun(int i, int n, double u0) { // calculates width of trapezoid and squares it
    return pow((i-1)*u0/(n-1),2);
 } 