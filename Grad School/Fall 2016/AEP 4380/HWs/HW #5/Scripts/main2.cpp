int main() {
	int istep, n = 16, i;
	double initvar[n] = {149.598e9, 228.0e9, 778.29e9/4.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 2.9786e4, 2.4127e4, 1.30588e4*sqrt(4.7), -3.0e1};
	double finalvar[16], t0 = 0.0, hold, h, nsecinday, totalsec, nsecinyear;
	void autosteprk4(double[], double[], int, double&, double, void(double[], double, double[], double));
	void nbody(double[], double, double[], double);
	ofstream fp;
    fp.open( "solarsystem2.dat" ); // open new file for output
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