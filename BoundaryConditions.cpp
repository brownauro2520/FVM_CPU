#include "BoundaryConditions.h"
#include "Solver.h" // Required for farField
#include <cmath>    // For std::sqrt, std::abs, std::isnan, std::isinf
#include <iostream> // For temporary error messages
#include <stdexcept>// For exceptions

BoundaryConditions::BoundaryConditions() {
    // Constructor
}

std::pair<std::vector<double>, double> BoundaryConditions::applyBoundaryCondition(
    const std::vector<double>& U_interior,
    const std::vector<double>& normal,
    int boundary_type_id,
    const std::vector<double>& u_far_field,
    Solver& solver
) {
    if (U_interior.size() != 4) {
        throw std::invalid_argument("U_interior must have 4 components.");
    }
    if (normal.size() != 2) {
        throw std::invalid_argument("Normal vector must have 2 components.");
    }

    // According to original residual.h:
    // if check == -1 -> farfield (use u_inf)
    // else -> airfoil boundary (use inviscid wall)
    // We'll map boundary_type_id == -1 to farField, and others to inviscidWall.
    // This might need to be more flexible later if more BC types are added.

    if (boundary_type_id == -1) { // Far-field condition
        if (u_far_field.size() != 4) {
             throw std::invalid_argument("u_far_field must have 4 components for far-field BC.");
        }
        return farField(U_interior, u_far_field, normal, solver);
    } else { // Inviscid wall condition (for any other type_id for now)
        return inviscidWall(U_interior, normal);
    }
}

// Implementation of inviscid wall condition (adapted from wall.h - inviscid function)
std::pair<std::vector<double>, double> BoundaryConditions::inviscidWall(
    const std::vector<double>& U_int, // Interior cell state U = [rho, rho*u, rho*v, E]
    const std::vector<double>& n_vec  // Normal vector n = [nx, ny]
) {
    // double gamma = 1.4; // Already a member gamma_
    double gmi = this->gamma_ - 1.0;

    double r = U_int[0];
    // Handle potential non-physical density from interior state
    if (r <= 1e-9) { // Threshold for too small or negative density
        // std::cerr << "Warning: Near-zero or negative density in inviscidWall: " << r << std::endl;
        // Option 1: Return zero flux / zero speed (problematic)
        // Option 2: Try to reconstruct a physical state (complex)
        // Option 3: Use a very small positive density
        r = 1e-9;
        // Or throw an error if this state is truly unrecoverable
        // throw std::runtime_error("Non-physical density encountered in inviscidWall.");
    }

    double ru = U_int[1];
    double rv = U_int[2];
    double E = U_int[3];

    double u = ru / r;
    double v = rv / r;

    // Velocity component normal to the wall
    double un = u * n_vec[0] + v * n_vec[1];

    // Pressure calculation from interior state: p = (gamma-1)*(E - 0.5*rho*(u^2+v^2))
    // For wall BC, the normal velocity component (un) is reflected.
    // The tangential component remains.
    // So, for the purpose of calculating pressure at the wall, we use the interior state's properties.
    // The original code's wall.h `inviscid` function:
    // q = sqrt(pow(U[1],2) + pow(U[2],2)) / r; -> magnitude of velocity |V|
    // ut = sqrt(pow(q,2) - pow(un,2)); -> magnitude of tangential velocity |Vt|
    // p = gmi*(U[3] - (0.5*U[0]*pow(ut,2))); -> This pressure calculation is unusual.
    // It seems to subtract tangential kinetic energy from total energy to get pressure.
    // Standard wall pressure is p_wall = p_interior.
    // Let's re-evaluate the pressure calculation.
    // p_interior = (gamma-1) * (E_interior - 0.5 * rho_interior * (u_interior^2 + v_interior^2))

    double p_int = gmi * (E - 0.5 * r * (u * u + v * v));
    if (p_int < 0) {
        // std::cerr << "Warning: Non-physical interior pressure in inviscidWall: " << p_int << std::endl;
        p_int = std::abs(p_int); // Or a small positive value
        if (p_int < 1e-9) p_int = 1e-9; // Ensure positive pressure for c_int calculation
    }
    double pressure_at_wall = p_int;


    // Flux vector F at the wall face.
    // F = [ rho*un_wall,
    //       rho*u*un_wall + p*nx,
    //       rho*v*un_wall + p*ny,
    //       (E+p)*un_wall ]
    // For an inviscid wall, normal velocity component un_wall = 0.
    std::vector<double> F_wall(4, 0.0);
    F_wall[0] = 0.0; // rho * 0
    F_wall[1] = pressure_at_wall * n_vec[0]; // p * nx
    F_wall[2] = pressure_at_wall * n_vec[1]; // p * ny
    F_wall[3] = 0.0; // (E+p) * 0

    // Maximum wave speed (signal speed) at the wall.
    // c = sqrt(gamma * p / rho)
    // smag = |un_interior| + c_interior
    // The original `wall.h` used `smag = abs(un) + c;` where un and c are from the interior cell.
    double c_int_sq = this->gamma_ * pressure_at_wall / r;
    double c_int = 0.0;
    if (c_int_sq < 0) {
        // This case should be rare if p_int and r are handled to be positive.
        // std::cerr << "Warning: Negative c_int_sq in inviscidWall: " << c_int_sq << std::endl;
        c_int_sq = std::abs(c_int_sq); // Take abs to prevent nan from sqrt, though this is non-physical
    }
    c_int = std::sqrt(c_int_sq);

    if (std::isnan(c_int) || std::isinf(c_int)) {
        // std::cerr << "Warning: Non-physical speed of sound c_int (NaN/Inf) in inviscidWall." << std::endl;
        // Fallback based on free stream, or a fixed small value. This is tricky.
        // For now, let's assume gamma, pressure_at_wall, and r are such that c_int is valid.
        // If p_at_wall is ~0 and r is ~0, c_int can be NaN.
        // If r is positive and p_at_wall is positive, c_int should be real.
        // Recheck handling of p_int and r.
        c_int = 1.0; // A default placeholder if calculation failed.
    }
    double smag = std::abs(un) + c_int;

    return std::make_pair(F_wall, smag);
}

// Implementation of far-field boundary condition
std::pair<std::vector<double>, double> BoundaryConditions::farField(
    const std::vector<double>& U_interior,
    const std::vector<double>& u_far_field_state,
    const std::vector<double>& normal,
    Solver& solver // Use the provided Riemann solver
) {
    // For far-field, the exterior state is u_far_field_state.
    // The solver calculates the flux between U_interior and u_far_field_state.
    return solver.calculateFlux(U_interior, u_far_field_state, normal);
}
