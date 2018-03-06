void chenlee(double var[], double t, double k[3], double h) {
	double x, y, z;
	x = var[0], y = var[1], z = var[2];
	k[0] = h*(5.0*x-y*z);
	k[1] = h*(x*z-10.0*y);
	k[2] = h*((1.0/3.0)*x*y-3.8*z);
	return;
}