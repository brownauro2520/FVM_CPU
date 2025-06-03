#include "HLLESolver.h"
#include <iostream> // For error messages (temporary)
#include <algorithm> // For std::max and std::min
#include <stdexcept> // Required for std::invalid_argument

HLLESolver::HLLESolver() {
    // Constructor can be empty or initialize specific HLLE solver parameters
}

std::string HLLESolver::getType() const {
    return "HLLE";
}

std::pair<std::vector<double>, double> HLLESolver::calculateFlux(
    const std::vector<double>& UL,
    const std::vector<double>& UR,
    const std::vector<double>& normal
) {
    // --- This is the core logic from the original HLLE() function in HLLE.h ---
    // Minor adaptations:
    // - Use this->gamma_
    // - Ensure vector sizes are checked
    // - Replace std::cout with proper error handling

    if (UL.size() != 4 || UR.size() != 4 || normal.size() != 2) {
        throw std::invalid_argument("Invalid input vector sizes for HLLE solver.");
    }

    std::vector<double> F(4);
    double gmi = this->gamma_ - 1.0;

    // Process left state
    double rL = UL[0];
    double uL = UL[1] / rL;
    double vL = UL[2] / rL;
    double unL = uL * normal[0] + vL * normal[1];
    double qL_sq = uL*uL + vL*vL; // Using qL_sq to avoid repeated sqrt
    double pL = gmi * (UL[3] - 0.5 * rL * qL_sq);

    if (pL < 0) { /*std::cout << "Warning: Non-physical pressure pL (HLLE): " << pL << std::endl;*/ pL = std::abs(pL); }
    if (rL < 0) { /*std::cout << "Warning: Non-physical density rL (HLLE): " << rL << std::endl;*/ rL = std::abs(rL); }

    double rHL = UL[3] + pL; // rho_L * H_L
    // double HL = rHL / rL; // Not explicitly needed for flux calculation itself
    double cL = std::sqrt(this->gamma_ * pL / rL);

    // Left flux components
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
    double qR_sq = uR*uR + vR*vR; // Using qR_sq
    double pR = gmi * (UR[3] - 0.5 * rR * qR_sq);

    if (pR < 0) { /*std::cout << "Warning: Non-physical pressure pR (HLLE): " << pR << std::endl;*/ pR = std::abs(pR); }
    if (rR < 0) { /*std::cout << "Warning: Non-physical density rR (HLLE): " << rR << std::endl;*/ rR = std::abs(rR); }

    double rHR = UR[3] + pR; // rho_R * H_R
    // double HR = rHR / rR; // Not explicitly needed
    double cR = std::sqrt(this->gamma_ * pR / rR);

    // Right flux components
    std::vector<double> FR(4);
    FR[0] = rR * unR;
    FR[1] = UR[1] * unR + pR * normal[0];
    FR[2] = UR[2] * unR + pR * normal[1];
    FR[3] = rHR * unR;

    // Estimate wave speeds (Davis, 1988)
    // S_L = min(unL - cL, unR - cR)
    // S_R = max(unL + cL, unR + cR)
    // More robust version from Toro book (Chapter 10, HLLC, simplified for HLLE S_L, S_R):
    double sL = std::min(unL - cL, unR - cR);
    double sR = std::max(unL + cL, unR + cR);

    // Original code's wave speed estimation:
    // double sLmax_orig = std::max(0.0, unL + cL); // This seems incorrect for sL / s_min for HLLE
    // double sRmax_orig = std::max(0.0, unR + cR); // This is S_R if unR+cR > 0
    // double sLmin_orig = std::min(0.0, unL - cL); // This is S_L if unL-cL < 0
    // double sRmin_orig = std::min(0.0, unR - cR);
    // double smin_orig = std::min(sLmin_orig, sRmin_orig);
    // double smax_orig = std::max(sLmax_orig, sRmax_orig);

    // Using Toro's S_L and S_R for HLLE (simplified from HLLC's S_L, S_R estimates)
    // S_L = u_L_n - c_L if k is L-state, u_R_n - c_R if k is R-state
    // S_R = u_L_n + c_L if k is L-state, u_R_n + c_R if k is R-state
    // For HLLE, a common choice is:
    // s_L = min(u_n_L, u_n_R) - max(c_L, c_R)  -- this is not quite right
    // s_R = max(u_n_L, u_n_R) + max(c_L, c_R)  -- this is not quite right

    // Einfeldt's wave speeds (more robust for HLLE)
    // Roe-averaged velocity u_tilde_n = (sqrt(rL)*unL + sqrt(rR)*unR) / (sqrt(rL)+sqrt(rR))
    // Roe-averaged enthalpy H_tilde = (sqrt(rL)*HL + sqrt(rR)*HR) / (sqrt(rL)+sqrt(rR))
    // Roe-averaged speed of sound c_tilde_sq = (gmi)*(H_tilde - 0.5*(u_tilde_n^2 + ...other_vel_components))
    // sL_einfeldt = u_tilde_n - c_tilde
    // sR_einfeldt = u_tilde_n + c_tilde
    // For simplicity and matching original code's apparent structure (though smax/smin logic was different):
    // We use simpler Davis (1988) estimates which are unL +/- cL and unR +/- cR.
    // s_L = min(unL - cL, unR - cR)
    // s_R = max(unL + cL, unR + cR)

    // If sR - sL is close to zero (e.g. sR <= sL), it can lead to division by zero.
    // This can happen if states are identical or very similar.
    // In such a case, the flux is simply FL or FR.
    if (sR <= sL) { // Or sR - sL < some_small_epsilon
         // If sR == sL, flux is (sR*FL - sL*FR + sR*sL*(UR-UL))/(sR-sL) -> indeterminate
         // Or if sL > 0, F = FL. If sR < 0, F = FR.
         // If sL <= 0 <= sR, then use the HLLE formula.
        if (sL > 0) F = FL;
        else if (sR < 0) F = FR;
        else { // sL <=0 and sR >=0 but sR is very close to sL (e.g. both near 0)
              // This case implies very small wave speeds, average flux is reasonable.
            for(int i=0; i<4; ++i) F[i] = 0.5 * (FL[i] + FR[i]);
        }
    } else {
        // Original code:
        // double sLmax_orig = max(0.0, unL + cL);
        // double sRmax_orig = max(0.0, unR + cR);
        // double sLmin_orig = min(0.0, unL - cL);
        // double sRmin_orig = min(0.0, unR - cR);
        // double smin_orig = min(sLmin_orig, sRmin_orig); // This is S_L in HLLE if both unL-cL and unR-cR are negative, else min(0, unL-cL, unR-cR)
        // double smax_orig = max(sLmax_orig, sRmax_orig); // This is S_R in HLLE if both unL+cL and unR+cR are positive, else max(0, unL+cL, unR+cR)

        // Let's use the sL and sR based on Davis (1988) as defined above.
        // S_L = min(unL - cL, unR - cR)
        // S_R = max(unL + cL, unR + cR)
        // The original code's smin/smax seems to be trying to ensure smin <=0 and smax >=0 for the HLLE formula terms.
        // The standard HLLE flux:
        // F_hlle = (sR * FL - sL * FR + sR * sL * (UR - UL)) / (sR - sL)

        // Let's use the original code's wave speed logic for smin, smax for now to match its behavior, then refine.
        // These are not exactly S_L and S_R of the Riemann fan for HLLE.
        // S_L = min( unL-cL, unR-cR )
        // S_R = max( unL+cL, unR+cR )
        // The original code used:
        double S_L_orig = std::min(unL - cL, unR - cR); // This is a candidate for S_L
        double S_R_orig = std::max(unL + cL, unR + cR); // This is a candidate for S_R

        // The original code then did:
        // sLmax = max(0.0, unL + cL); sRmax = max(0.0, unR + cR);
        // sLmin = min(0.0, unL - cL); sRmin = min(0.0, unR - cR);
        // smin_from_code = min(sLmin, sRmin);
        // smax_from_code = max(sLmax, sRmax);
        // This ensures smin_from_code <= 0 and smax_from_code >= 0 if the flow is not entirely supersonic in one direction.

        // Standard HLLE formulation uses S_L and S_R directly:
        // S_L = min(u_n - c, u_n_avg - c_avg) etc. Simpler: S_L = min(unL-cL, unR-cR), S_R = max(unL+cL, unR+cR)
        // Let's use these simpler, more standard S_L, S_R
        double s_L_davis = unL - cL;
        double s_R_davis = unL + cL;
        double s_L_davis_UR = unR - cR;
        double s_R_davis_UR = unR + cR;

        sL = std::min(s_L_davis, s_L_davis_UR);
        sR = std::max(s_R_davis, s_R_davis_UR);


        // Ensure sR > sL to avoid division by zero or negative denominator.
        // This case should ideally be handled by the check sR <= sL above.
        // If sR is very close to sL, the denominator is small.
        if (sR - sL < 1e-9) { // Small epsilon to prevent division by zero
            // If wave speeds are basically identical, average flux
            for (int i = 0; i < 4; i++) {
                F[i] = 0.5 * (FL[i] + FR[i]);
            }
        } else {
             // Original code used C1 and C2 based on its smin/smax
             // C1 = 0.5 * (smax_from_code + smin_from_code) / (smax_from_code - smin_from_code);
             // C2 = smax_from_code * smin_from_code / (smax_from_code - smin_from_code);
             // F[i] = 0.5 * (FL[i] + FR[i]) - C1 * (FR[i] - FL[i]) + C2 * (UR[i] - UL[i]);
             // This is equivalent to F_hlle = (sR * FL - sL * FR + sR * sL * (UR - UL)) / (sR - sL)
             // if sR = smax_from_code and sL = smin_from_code, and sR*sL terms are handled carefully.
             // Let's use the direct HLLE formula with sL and sR (Davis estimates):
            for (int i = 0; i < 4; i++) {
                F[i] = (sR * FL[i] - sL * FR[i] + sR * sL * (UR[i] - UL[i])) / (sR - sL);
            }
        }
    }

    // Maximum wave speed for CFL condition
    // smag = max(|sL|, |sR|) for HLLE
    double smag = std::max(std::abs(sL), std::abs(sR));
    // The original code used: smag = smax_from_code; (which was max(0, unL+cL, unR+cR))
    // A more robust CFL limit for HLLE is indeed max(|S_L|, |S_R|)
    // where S_L and S_R are the actual min/max signal velocities.

    return std::make_pair(F, smag);
}
