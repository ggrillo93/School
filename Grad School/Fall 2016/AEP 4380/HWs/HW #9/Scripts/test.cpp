#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>
#include "arrayt.hpp"
#include "nr3.h"
#include "ran.h"
#define ARRAYT_BOUNDS_CHECK

using namespace std;

int main() {
    int naa = 45, i, t;
    struct Ran myrand(17);
    void aamodifier(arrayt<int>&, int&, int);
    double sep(int&, int&, int&, int&);
    arrayt<int> aalist(naa, 3);
    for (i = 0; i < naa; i++) { // initialize amino acids horizontally
        aalist(i,0) = i; // initial x-coordinate
        aalist(i,1) = 0; // initial y-coordinate
        aalist(i,2) = (myrand.int64() % 20); // type of amino acid
    }
    sep(aalist(0,0), aalist(naa-1, 0), aalist(0, 1), aalist(naa-1, 1));
    return(EXIT_SUCCESS);
}

double energy(arrayt<int>& aalist, arrayt<double>& enmat) {
    int naa = aalist.n1(), i, j, x1, x2, y1, y2, type1, type2;
    double energy = 0.0;
    int sep(int&, int&, int&, int&);
    for (i = 0; i < naa; i++) for (j = 0; j < i; j++) {
        if (j != i+1 && j != i-1) { // no need to check covalent bonds
            x1 = aalist(i, 0);
            x2 = aalist(j, 0);
            y1 = aalist(i, 1);
            y2 = aalist(j, 1);
            if (sep(x1, x2, y1, y2) == 1) {
                type1 = aalist(i, 2);
                type2 = aalist(j, 2);
                energy += enmat(type1, type2);
            }
        }
    }
    return energy;
}

int sep(int& x1, int& x2, int& y1, int& y2) {
    cout << (abs(x1 - x2) + abs(y1 - y2)) << endl;
    return (abs(x1 - x2) + abs(y1 - y2));
}

void aamodifier(arrayt<int>& aalist, int& naa, int rn2) {
    bool inarray(arrayt<int>&, int&, int&), allowed(int&, int&, int&, int&, int&, int&);
    void randomize(arrayt<int>&, int&);
    int p1x, p1y, p2x, p2y, newx, newy, oldx, oldy, i, j, aan, rn;
    arrayt<int> ranaa(naa), ranmov(4);
    for (i = 0; i < naa; i++) { ranaa(i) = i; }
    for (i = 0; i < 4; i++) { ranmov(i) = i; }
    randomize(ranaa, rn2);
    randomize(ranmov, rn2);
    for (i = 0; i < naa; i++) for (j = 0; j < 8; j++) {
        aan = ranaa(i);
        rn = ranmov(j);
        oldx = aalist(aan, 0);
        oldy = aalist(aan, 1);
        if (rn == 0) {
            newx = oldx - 1;
            newy = oldy + 1;
            if (aan == 0) { // special case 1
                p1x = aalist(aan + 1, 0);
                p1y = aalist(aan + 1, 1);
                if (inarray(aalist, newx, newy) == false && sep(newx, p1x, newy, p1y) == 1) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
            if (aan == naa-1) { // special case 2
                p1x = aalist(aan - 1, 0);
                p1y = aalist(aan - 1, 1);
                if (inarray(aalist, newx, newy) == false && sep(newx, p1x, newy, p1y) == 1) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
            else {
                p1x = aalist(aan + 1, 0);
                p1y = aalist(aan + 1, 1);
                p2x = aalist(aan - 1, 0);
                p2y = aalist(aan - 1, 1);
                if (inarray(aalist, newx, newy) == false && allowed(newx, newy, p1x, p1y, p2x, p2y)) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
        }
        if (rn == 1) {
            newx = oldx + 1;
            newy = oldy + 1;
            if (aan == 0) { // special case 1
                p1x = aalist(aan + 1, 0);
                p1y = aalist(aan + 1, 1);
                if (inarray(aalist, newx, newy) == false && sep(newx, p1x, newy, p1y) == 1) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
            if (aan == naa-1) { // special case 2
                p1x = aalist(aan - 1, 0);
                p1y = aalist(aan - 1, 1);
                if (inarray(aalist, newx, newy) == false && sep(newx, p1x, newy, p1y) == 1) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
            else {
                p1x = aalist(aan + 1, 0);
                p1y = aalist(aan + 1, 1);
                p2x = aalist(aan - 1, 0);
                p2y = aalist(aan - 1, 1);
                if (inarray(aalist, newx, newy) == false && allowed(newx, newy, p1x, p1y, p2x, p2y)) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
        }
        if (rn == 2) {
            newx = oldx + 1;
            newy = oldy - 1;
            if (aan == 0) { // special case 1
                p1x = aalist(aan + 1, 0);
                p1y = aalist(aan + 1, 1);
                if (inarray(aalist, newx, newy) == false && sep(newx, p1x, newy, p1y) == 1) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
            if (aan == naa-1) { // special case 2
                p1x = aalist(aan - 1, 0);
                p1y = aalist(aan - 1, 1);
                if (inarray(aalist, newx, newy) == false && sep(newx, p1x, newy, p1y) == 1) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
            else {
                p1x = aalist(aan + 1, 0);
                p1y = aalist(aan + 1, 1);
                p2x = aalist(aan - 1, 0);
                p2y = aalist(aan - 1, 1);
                if (inarray(aalist, newx, newy) == false && allowed(newx, newy, p1x, p1y, p2x, p2y)) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
        }
        if (rn == 3) {
            newx = oldx - 1;
            newy = oldy - 1;
            if (aan == 0) { // special case 1
                p1x = aalist(aan + 1, 0);
                p1y = aalist(aan + 1, 1);
                if (inarray(aalist, newx, newy) == false && sep(newx, p1x, newy, p1y) == 1) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
            if (aan == naa-1) { // special case 2
                p1x = aalist(aan - 1, 0);
                p1y = aalist(aan - 1, 1);
                if (inarray(aalist, newx, newy) == false && sep(newx, p1x, newy, p1y) == 1) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
            else {
                p1x = aalist(aan + 1, 0);
                p1y = aalist(aan + 1, 1);
                p2x = aalist(aan - 1, 0);
                p2y = aalist(aan - 1, 1);
                if (inarray(aalist, newx, newy) == false && allowed(newx, newy, p1x, p1y, p2x, p2y)) {
                    aalist(aan, 0) = newx;
                    aalist(aan, 1) = newy;
                    goto end_loop;
                }
            }
        }
    }
    end_loop:
    return;
}

bool inarray(arrayt<int>& aalist, int& newx, int& newy) {
    int size = aalist.n1(), i, x, y;
    bool r = false;
    for (i = 0; i < size; i++) {
        x = aalist(i, 0);
        y = aalist(i, 1);
        if (newx == x && newy == y) { r = true; }
    }
    return r;
}

bool allowed(int& newx, int& newy, int& p1x, int& p1y, int& p2x, int& p2y) {
    int sep(int&, int&, int&, int&);
    int sep1, sep2;
    sep1 = sep(newx, p1x, newy, p1y);
    sep2 = sep(newx, p2x, newy, p2y);
    if (sep1 > 1 || sep2 > 1) { return false; }
    else { return true; }
}

void randomize(arrayt<int>& list, int& rn) {
    int size, i, rn1, rn2;
    struct Ran myrand2(rn);
    size = list.n();
    for (i = 0; i < size; i++) {
        rn1 = myrand2.int64() % size;
        rn2 = myrand2.int64() % size;
        SWAP(list(rn1), list(rn2));
    }
    return;
}
