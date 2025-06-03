#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>
#include <cstdlib> // For std::atof

#include "wall.h"
#include "roe.h"
#include "HLLE.h"
#include "extrafunctions.h"
#include "residual.h"

using namespace std;

int main(int argc, char* argv[]) {
    double M = 0.25; //subsonic = 0.25
    double CFL = 1.2; // Default CFL value

    // Check if a CFL value was provided as a system argument
    if (argc > 1) {
        CFL = atof(argv[1]); // Convert the first argument to a double and set as CFL
    }

    double alpha = 8 * M_PI / 180;
    double gamma = 1.4;
    double gmi = gamma - 1;
    vector<double> u_inf = { 1, M * cos(alpha), M * sin(alpha), (1 / (gmi * gamma)) + (pow(M,2) / 2) };

    vector<vector<double>> E2N = readmatrix("matrices/E2N.txt");
    vector<vector<double>> I2E = readmatrix("matrices/I2E.txt");
    vector<vector<double>> B2E = readmatrix("matrices/B2E.txt");
    vector<vector<double>> In = readmatrix("matrices/In.txt");
    vector<vector<double>> Bn = readmatrix("matrices/Bn.txt");
    vector<double> Il = readvector("matrices/Il.txt");
    vector<double> Bl = readvector("matrices/Bl.txt");
    vector<double> A = readvector("matrices/A.txt");
    int N = E2N.size();

    vector<double> U(N * 4);
    for (int i = 0; i < N * 4; i += 4) {
        U[i] = u_inf[0]; U[i + 1] = u_inf[1]; U[i + 2] = u_inf[2]; U[i + 3] = u_inf[3];
    }

    // Time stepping variables
    vector<double> store_res;
    double tolerance = pow(10, -5);

    // Allocate memory for k stages and intermediate variables
    vector<double> k1(N * 4), k2(N * 4), k3(N * 4), k4(N * 4);
    vector<double> u_temp(N * 4);
    double total_res;
    vector<double> dT(N);
    vector<double> R(N * 4);
    auto start = chrono::high_resolution_clock::now();

    for (int iter = 0; iter < 1000000; iter++) {
        fill(dT.begin(), dT.end(), 0.0);
        fill(R.begin(), R.end(), 0.0);
        total_res = 0;

        // Stage 1: Calculate k1
        residual(R, dT, U, I2E, B2E, In, Bn, Il, Bl, A, u_inf, N, CFL);
        for (int i = 0; i < N * 4; i++) {
            k1[i] = R[i] / (-A[floor(i / 4)]);
        }

        // Stage 2: Calculate k2
        for (int i = 0; i < N * 4; i++) {
            u_temp[i] = U[i] + 0.5 * k1[i] * dT[floor(i / 4)];
        }
        residual(R, dT, u_temp, I2E, B2E, In, Bn, Il, Bl, A, u_inf, N, CFL);
        for (int i = 0; i < N * 4; i++) {
            k2[i] = R[i] / (-A[floor(i / 4)]);
        }

        // Stage 3: Calculate k3
        for (int i = 0; i < N * 4; i++) {
            u_temp[i] = U[i] + 0.5 * k2[i] * dT[floor(i / 4)];
        }
        residual(R, dT, u_temp, I2E, B2E, In, Bn, Il, Bl, A, u_inf, N, CFL);
        for (int i = 0; i < N * 4; i++) {
            k3[i] = R[i] / (-A[floor(i / 4)]);
        }

        // Stage 4: Calculate k4
        for (int i = 0; i < N * 4; i++) {
            u_temp[i] = U[i] + k3[i] * dT[floor(i / 4)];
        }
        residual(R, dT, u_temp, I2E, B2E, In, Bn, Il, Bl, A, u_inf, N, CFL);
        for (int i = 0; i < N * 4; i++) {
            k4[i] = R[i] / (-A[floor(i / 4)]);
        }

        // Update solution
        for (int i = 0; i < N * 4; i++) {
            U[i] += (1.0 / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * dT[floor(i / 4)];
        }

        // Compute total residual
        for (int i = 0; i < N * 4; i++) {
            total_res += abs(R[i]);
        }

        cout << "Residual " << total_res << " iter " << iter << endl;

        // Storing residual
        store_res.push_back(abs(total_res));

        // Check for convergence
        if (abs(total_res) < tolerance) {
            savefilev(U, "sols/U.txt");
            auto stop = chrono::high_resolution_clock::now();
            auto time = chrono::duration_cast<chrono::milliseconds>(stop - start);
            double dt = time.count();
            dt = dt / 1000;
            cout << "Time to Reach Solution: " << dt << " seconds \n";
            break;
        }
    }

    return 0;
}
