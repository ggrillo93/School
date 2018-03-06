/* AEP 4380 HW#9a
  Monte Carlo Calculations
  Run on core i7 with g++ 5.4.0 (Ubuntu)
  Gianfranco Grillo 11/20/2016
*/

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>
#include "arrayt.hpp"
#include "nr3.h"
#include "ran.h"

using namespace std;

int main() {
    int npts, npairs, nbins, n = 0, i;
    double rn, sqrsum;
    ofstream fp, sp;
    cout << "Enter number of points: ";
    cin >> npts;
    cout << "Enter number of bins: ";
    cin >> nbins;
    cout << "Enter number of pairs: ";
    cin >> npairs;
    struct Ran myrand(17);
    arrayt<int> hist(nbins);
    arrayt<double> pairs(npairs,2), histranges(nbins+1);
    fp.open("histogram.dat");
    sp.open("pairs.dat");
    for (i = 0.0; i < nbins+1.0; i++) { histranges(i) = i; }
    for (i = 0; i < nbins; i++) { hist(i) = 0; }
    while (n < npts) {
        rn = myrand.doub();
        for (i = 0; i < histranges.n(); i++) {
            if (rn > (histranges(i)/nbins) && rn <= (histranges(i+1)/nbins)) {
                hist(i)++;
                break;
            }
        }
        n++;
    }
    sqrsum = 0.0;
    for (i = 0; i < npairs; i = i+2) {
        pairs(i, 0) = myrand.doub();
        pairs(i, 1) = myrand.doub();
        sp << pairs(i, 0) << setw(15) << pairs(i, 1) << endl;
    }
    for (i = 0; i < nbins; i++) {
        fp << hist(i) << endl;
    }
    fp.close();
    sp.close();
    return(EXIT_SUCCESS);
}
