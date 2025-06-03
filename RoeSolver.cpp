#include "RoeSolver.h"
#include <iostream> // For error messages (temporary)
#include <numeric>  // For std::iota or other utilities if needed
#include <stdexcept> // Required for std::invalid_argument

RoeSolver::RoeSolver() {
    // Constructor can be empty or initialize specific Roe solver parameters
}

std::string RoeSolver::getType() const {
    return "Roe";
}

// Helper function implementations (can be made static if they don't access member variables)
std::vector<double> RoeSolver::absolute_value(const std::vector<double>& v) {
    std::vector<double> abs_v(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        abs_v[i] = std::abs(v[i]);
    }
    return abs_v;
}

std::vector<double> RoeSolver::subtract(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        // Handle error: vectors must have the same size
        throw std::invalid_argument("Vectors must have the same size for subtraction.");
    }
    std::vector<double> result(v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

std::pair<std::vector<double>, double> RoeSolver::calculateFlux(
    const std::vector<double>& UL,
    const std::vector<double>& UR,
    const std::vector<double>& normal
) {
    // --- This is the core logic from the original roe() function in roe.h ---
    // Minor adaptations:
    // - Use this->gamma_ or just gamma_
    // - Ensure vector sizes are checked or handled robustly
    // - Replace std::cout with proper error handling or logging if necessary

    if (UL.size() != 4 || UR.size() != 4 || normal.size() != 2) {
        throw std::invalid_argument("Invalid input vector sizes for Roe solver.");
    }

    double gmi = this->gamma_ - 1.0;

    // Process left state
    double rL = UL[0];
    double uL = UL[1] / rL;
    double vL = UL[2] / rL;
    double unL = uL * normal[0] + vL * normal[1];
    double qL_sq = uL * uL + vL * vL; // qL^2
    double pL = gmi * (UL[3] - 0.5 * rL * qL_sq);

    if (pL < 0) { /*std::cout << "Warning: Non-physical pressure pL: " << pL << std::endl;*/ pL = std::abs(pL); }
    if (rL < 0) { /*std::cout << "Warning: Non-physical density rL: " << rL << std::endl;*/ rL = std::abs(rL); }


    double rHL = UL[3] + pL; // rho_L * H_L
    double HL = rHL / rL;
    double cL = std::sqrt(this->gamma_ * pL / rL);

    // Left flux
    std::vector<double> FL(4);
    FL[0] = rL * unL;
    FL[1] = UL[1] * unL + pL * normal[0];
    FL[2] = UL[2] * unL + pL * normal[1];
    FL[3] = rHL * unL;

    // Process right state
    double rR = UR[0];
    double uR = UR[1] / rR;
    double vR = UR[2] / rR;
    double unR = uR * normal[0] + vR * normal[1];
    double qR_sq = uR*uR + vR*vR; // qR^2
    double pR = gmi * (UR[3] - 0.5 * rR * qR_sq);

    if (pR < 0) { /*std::cout << "Warning: Non-physical pressure pR: " << pR << std::endl;*/ pR = std::abs(pR); }
    if (rR < 0) { /*std::cout << "Warning: Non-physical density rR: " << rR << std::endl;*/ rR = std::abs(rR); }

    double rHR = UR[3] + pR; // rho_R * H_R
    double HR = rHR / rR;
    double cR = std::sqrt(this->gamma_ * pR / rR);

    // Right flux
    std::vector<double> FR(4);
    FR[0] = rR * unR;
    FR[1] = UR[1] * unR + pR * normal[0];
    FR[2] = UR[2] * unR + pR * normal[1];
    FR[3] = rHR * unR;

    // Difference in states
    std::vector<double> du = subtract(UR, UL);

    // Roe average
    double di = std::sqrt(rR / rL); // sqrt(rho_R / rho_L)
    double d1 = 1.0 / (1.0 + di);
    double ui = (di * uR + uL) * d1;
    double vi = (di * vR + vL) * d1;
    double Hi = (di * HR + HL) * d1;
    double af = 0.5 * (ui * ui + vi * vi); // average u_tilde^2 / 2
    double ucp = ui * normal[0] + vi * normal[1]; // u_tilde_n
    double c2 = gmi * (Hi - af); // c_tilde^2

    if (c2 < 0) {
        //std::cout << "Warning: Non-physical c2 (Roe): " << c2 << std::endl;
        c2 = std::abs(c2); // Or handle more robustly, e.g. use average of cL, cR
    }
    double ci = std::sqrt(c2); // c_tilde
    double ci1 = 1.0 / ci;

    // Eigenvalues
    std::vector<double> l(3);
    l[0] = ucp + ci;
    l[1] = ucp - ci;
    l[2] = ucp;

    // Entropy fix (Harten's fix)
    double epsilon = ci * 0.1; // Typical epsilon value
    for (int i = 0; i < 3; i++) {
        if (std::abs(l[i]) < epsilon) {
            l[i] = 0.5 * (l[i] * l[i] / epsilon + epsilon);
        }
    }

    std::vector<double> l_abs = absolute_value(l); // Store absolute values
    double l3_abs = l_abs[2];


    // Average and half-difference of 1st and 2nd eigs (absolute values)
    double s1 = 0.5 * (l_abs[0] + l_abs[1]);
    double s2 = 0.5 * (l_abs[0] - l_abs[1]);

    // Left eigenvector product generators
    double G1 = gmi * (af * du[0] - ui * du[1] - vi * du[2] + du[3]);
    double G2 = -ucp * du[0] + du[1] * normal[0] + du[2] * normal[1];

    // Required functions of G1 and G2
    double C1 = G1 * (s1 - l3_abs) * ci1 * ci1 + G2 * s2 * ci1;
    double C2 = G1 * s2 * ci1 + G2 * (s1 - l3_abs);

    // Flux assembly
    std::vector<double> F(4);
    F[0] = 0.5 * (FL[0] + FR[0]) - 0.5 * (l3_abs * du[0] + C1);
    F[1] = 0.5 * (FL[1] + FR[1]) - 0.5 * (l3_abs * du[1] + C1 * ui + C2 * normal[0]);
    F[2] = 0.5 * (FL[2] + FR[2]) - 0.5 * (l3_abs * du[2] + C1 * vi + C2 * normal[1]);
    F[3] = 0.5 * (FL[3] + FR[3]) - 0.5 * (l3_abs * du[3] + C1 * Hi + C2 * ucp);

    // Max wave speed
    double smag = 0.0;
    // The original code had:
    // auto max_s = max_element(l.begin(), l.end()); // This was looking at non-absolute eigenvalues for smag
    // double smag = *max_s;
    // Corrected to use absolute eigenvalues for stability (max of |ucp|+c, |ucp|-c, |ucp|)
    // Or more simply, max of |ucp|+ci based on Roe averages, or max of |unL|+cL and |unR|+cR
    smag = std::max(std::abs(unL) + cL, std::abs(unR) + cR);
    // A common definition for Roe scheme's max wave speed for CFL is max(|ucp|+ci, |ucp-ci|, |ucp|)
    // smag = std::max({std::abs(l[0]), std::abs(l[1]), std::abs(l[2])});
    // Let's use the max of physical wave speeds at L and R states as a robust measure for CFL
    smag = std::max(std::abs(unL) + cL, std::abs(unR) + cR);


    return std::make_pair(F, smag);
}
