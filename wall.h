#ifndef WALL_H_INCLUDED
#define WALL_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

pair<vector<double>, double> inviscid(vector<double> &U, vector<double> &n) {
double gamma = 1.4;
double gmi = gamma-1;
double r = U[0];
double u = U[1]/U[0];
double v = U[2]/U[0];
double un = u*n[0] + v*n[1];
double q = sqrt(pow(U[1],2) + pow(U[2],2)) / r;
double ut = sqrt(pow(q,2) - pow(un,2));
double p = gmi*(U[3] - (0.5*U[0]*pow(ut,2)));
if (p < 0 || r < 0){
    cout << "Non-physical State! Inviscid";
}
if (p < 0){p = -p;}
if (r < 0){r = -r;}
double rH = U[3] + p;
double c = sqrt(gamma*p / r);
double smag = abs(un) + c;

vector<double> F(4,0);
F[1] = p*n[0];
F[2] = p*n[1];

return make_pair(F, smag);
}

#endif // WALL_H_INCLUDED
