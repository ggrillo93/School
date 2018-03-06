void matmult(arrayt<double>& a, arrayt<double>& b, arrayt<double>& c) {
	int i, j, k;
    if (a.n2() != b.n1()) {
        cout << "Matrix multiplication can not be performed." << endl;
        return;
    }
    for (i = 0; i < a.n1(); i++) {
        for (j = 0; j < b.n2(); j++) {
            double temp(0.0);
            for (k = 0; k < a.n2(); k++) {
                temp += a(i,k)*b(k,j);
            }
            c(i,j) = temp;
        }
    }
    return;
}
