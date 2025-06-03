#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace std;


vector<double> absolute_value(vector<double> v) {
  vector<double> abs_v(v.size());
  for (int i = 0; i < v.size(); i++) {
    abs_v[i] = abs(v[i]);
  }
  return abs_v;
}

vector<double> subtract(const vector<double> &v1, const vector<double> &v2) {
  vector<double> result(v1.size());
  for (int i = 0; i < v1.size(); i++) {
    result[i] = v1[i] - v2[i];
  }
  return result;
}

pair<vector<double>, double> roe(vector<double> &UL, vector<double> &UR, vector<double> &n, double gamma) {

  double gmi = gamma - 1;

  //process left state
  double rL = UL[0];
  double uL = UL[1]/rL;
  double vL = UL[2]/rL;
  double unL = uL*n[0] + vL*n[1];
  double qL = sqrt(pow(UL[1],2) + pow(UL[2],2))/rL;
  double pL = (gamma-1)*(UL[3] - 0.5*rL*pow(qL,2));

  //if (pL < 0 || rL < 0) {cout << "Non-physical state!1" << endl;}
  if (pL < 0){pL = -pL;}
  if (rL < 0){rL = -rL;}

  double rHL = UL[3] + pL;
  double HL = rHL/rL;
  double cL = sqrt(gamma*pL/rL);

  //left flux
  vector<double> FL(4);
  FL[0] = rL*unL;
  FL[1] = UL[1]*unL + pL*n[0];
  FL[2] = UL[2]*unL + pL*n[1];
  FL[3] = rHL*unL;

  //process right state
  double rR = UR[0];
  double uR = UR[1]/rR;
  double vR = UR[2]/rR;
  double unR = uR*n[0] + vR*n[1];
  double qR = sqrt(pow(UR[1],2) + pow(UR[2],2))/rR;
  double pR = (gamma-1)*(UR[3] - 0.5*rR*pow(qR,2));

 // if ((pR < 0) || (rR < 0)) { cout << "Non-physical state!2" << endl; }
  if (pR < 0){pR = -pR;}
  if (rR < 0){rR = -rR;}

  double rHR = UR[3] + pR;
  double HR = rHR/rR;
  double cR = sqrt(gamma*pR/rR);

  //right flux
  vector<double> FR(4);
  FR[0] = rR*unR;
  FR[1] = UR[1]*unR + pR*n[0];
  FR[2] = UR[2]*unR + pR*n[1];
  FR[3] = rHR*unR;


  //difference in states
  //vector<double> du(UL.size());
  vector<double> du = subtract(UR,UL);

  //Roe average
  double di = sqrt(rR/rL);
  double d1 = 1.0/(1.0+di);
  double ui = (di*uR + uL)*d1;
  double vi = (di*vR + vL)*d1;
  double Hi = (di*HR + HL)*d1;
  double af = 0.5*(ui*ui+vi*vi);
  double ucp = ui*n[0] + vi*n[1];
  double c2 = gmi*(Hi - af);
  if (c2 < 0){cout << "Non-physical state!3"; c2 = -c2;}
  double ci = sqrt(c2);
  double ci1 = 1.0/ci;

  //eigenvalues
  vector<double> l(3);
  l[0] = ucp+ci;
  l[1] = ucp-ci;
  l[2] = ucp;

  //entropy fix
  double epsilon = ci*0.1;
  for (int i = 0; i < 3; i++){
    if ((l[i]< epsilon) & (l[i]> -epsilon)){
        l[i] = 0.5*(epsilon + l[i]*l[i]/epsilon);
    }
  }

  l = absolute_value(l);
  double l3 = l[2];

  //average and half-difference of 1st and 2nd eigs
  double s1 = 0.5*(l[0] + l[1]);
  double s2 = 0.5*(l[0] - l[1]);

  //left eigenvector product generators
  double G1 = gmi*(af*du[0] - ui*du[1] - vi*du[2] + du[3]);
  double G2 = -ucp*du[0]+du[1]*n[0]+du[2]*n[1];

  //required functions of G1 and G2
  double C1 = G1*(s1-l3)*ci1*ci1 + G2*s2*ci1;
  double C2 = G1*s2*ci1 + G2*(s1-l3);

  //flux assembly
  vector<double> F(4);
  F[0] = 0.5*(FL[0]+FR[0])-0.5*(l3*du[0] + C1);
  F[1] = 0.5*(FL[1]+FR[1])-0.5*(l3*du[1] + C1*ui + C2*n[0]);
  F[2] = 0.5*(FL[2]+FR[2])-0.5*(l3*du[2] + C1*vi + C2*n[1]);
  F[3] = 0.5*(FL[3]+FR[3])-0.5*(l3*du[3] + C1*Hi + C2*ucp);

 // max wave speed
  auto max_s = max_element(l.begin(), l.end());
  double smag = *max_s;

  return make_pair(F, smag);
}