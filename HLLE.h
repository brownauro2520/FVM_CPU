#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <cmath>

using namespace std;

pair<vector<double>, double> HLLE(vector<double>& UL, vector<double>& UR, vector<double>& n, double gamma) {
    vector<double> F = { 0,0,0,0 };

    // Process left state
    double rL = UL[0];
    double uL = UL[1] / rL;
    double vL = UL[2] / rL;
    double unL = uL * n[0] + vL * n[1];
    double qL = std::sqrt(pow(UL[1], 2) + pow(UL[2], 2)) / rL;
    double pL = (gamma - 1) * (UL[3] - 0.5 * rL * pow(qL, 2));

    if (pL < 0 || rL < 0) {std::cout << "Non-physical state!" << std::endl;}
    if (pL < 0) { pL = -pL; }
    if (rL < 0) { rL = -rL; }

    double rHL = UL[3] + pL;
    double HL = rHL / rL;
    double cL = std::sqrt(gamma * pL / rL);

    //left flux
    std::vector<double> FL = { 0,0,0,0 };
    FL[0] = rL * unL;
    FL[1] = UL[1] * unL + pL * n[0];
    FL[2] = UL[2] * unL + pL * n[1];
    FL[3] = rHL * unL;

    //process right state
    double rR = UR[0];
    double uR = UR[1] / rR;
    double vR = UR[2] / rR;
    double unR = uR * n[0] + vR * n[1];
    double qR = sqrt(pow(UR[1], 2) + pow(UR[2], 2)) / rR;
    double pR = (gamma - 1) * (UR[3] - 0.5 * rR * pow(qR, 2));

    if ((pR < 0) || (rR < 0)) { cout << "Non-physical state!" << endl; }
    if (pR < 0) { pR = -pR; }
    if (rR < 0) { rR = -rR; }

    double rHR = UR[3] + pR;
    double HR = rHR / rR;
    double cR = sqrt(gamma * pR / rR);

    // right flux
    std::vector<double> FR = { 0,0,0,0 };
    FR[0] = rR * unR;
    FR[1] = UR[1] * unR + pR * n[0];
    FR[2] = UR[2] * unR + pR * n[1];
    FR[3] = rHR * unR;


    double sLmax = max(0.0, unL + cL);
    double sRmax = max(0.0, unR + cR);
    double sLmin = min(0.0, unL - cL);
    double sRmin = min(0.0, unR - cR);
    double smin = min(sLmin, sRmin);
    double smax = max(sLmax, sRmax);
    double C1 = 0.5 * (smax + smin) / (smax - smin);
    double C2 = smax * smin / (smax - smin);


    for (int i = 0; i < 4; i++) 
    {
        F[i] = 0.5 * (FL[i] + FR[i]) - C1 * (FR[i] - FL[i]) + C2 * (UR[i] - UL[i]);
    }
    double smag = smax;
    return make_pair(F,smag);

}