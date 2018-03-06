int main() {
    int i, npts = 607, npar = 9, j, k, t;
    double pi = 4.0 * atan(1.0), Flk, bl, chisqr;
    arrayt<int> tlist(npts);
    arrayt<double> co2list(npts), blist(npar), fit(npts), sqrerrlist(npts), reslist(npts);
    arrayt<double> Fmatrix(npar, npar), funcm(npts,npar);
    void readfile(arrayt<int>&, arrayt<double>&); // create t, co2, arrays from file
    void gaussj(arrayt<double>&, arrayt<double>&);
    ofstream fp, sp;
    fp.open("output7.dat");
    sp.open("output8.dat");
    readfile(tlist, co2list);
    for (i = 0; i < npts; i++) { // initialize function arrays
        t = tlist(i);
        funcm(i,0) = sin(pi*t/6.0);
        funcm(i,1) = cos(pi*t/6.0);
        funcm(i,2) = sin(pi*t/3.0);
        funcm(i,3) = cos(pi*t/3.0);
        funcm(i,4) = t*t;
        funcm(i,5) = t; // don't really need this but makes life easier
        funcm(i,6) = 1.0;
        funcm(i,7) = cos(pi*t/150);
        funcm(i,8) = sin(pi*t/150);
        sqrerrlist(i) = (0.002*co2list(i))*(0.002*co2list(i)); // create error array
    }
    for (i = 0; i < npar; i++) for (j = 0; j < i+1; j++) { // only need to calculate half of matrix
        Flk = 0;
        for (k = 0; k < npts; k++) {
            Flk = Flk + funcm(k,i)*funcm(k,j)/sqrerrlist(k);
        }
        Fmatrix(i,j) = Flk;
        if (i != j) {
            Fmatrix(j,i) = Flk;
        }
    }
    for (i = 0; i < npar; i++) { // generate b vector
        bl = 0;
        for (k = 0; k < npts; k++) {
            bl = bl + co2list(k)*funcm(k,i)/sqrerrlist(k);
        }
        blist(i) = bl;
    }
    gaussj(Fmatrix, blist);
    chisqr = 0.0;
    for (i = 0; i < npts; i++) {
        fit(i) = 0;
        for (j = 0; j < npar; j++) {
            fit(i) = fit(i) + blist(j)*funcm(i, j);
            }
        reslist(i) = co2list(i)-fit(i);
        chisqr = chisqr + reslist(i)*reslist(i)/sqrerrlist(i);
        fp << tlist(i) << setw(15) << co2list(i) << setw(15) << fit(i) << setw(15) << reslist(i) << endl;
    }
    chisqr = chisqr/(npts-npar);
    cout << chisqr << endl;
    for (i = 0; i < npar; i++) {
        sp << Fmatrix(i,i) << setw(15) << blist(i) << endl;
    }
    fp.close();
    return (EXIT_SUCCESS);
}
