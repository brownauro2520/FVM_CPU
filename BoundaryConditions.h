#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <vector>
#include <string>
#include <utility> // For std::pair

// Forward declarations
class Solver; // To use a Riemann solver for some BCs like far-field

class BoundaryConditions {
public:
    BoundaryConditions();
    ~BoundaryConditions() = default;

    // Applies the appropriate boundary condition and returns flux and max wave speed.
    //
    // Parameters:
    //   U_interior: State vector of the interior cell adjacent to the boundary.
    //   normal: Normal vector of the boundary face (pointing outwards from the domain).
    //   boundary_type_id: An identifier for the type of boundary
    //                     (e.g., -1 for far-field, other positive integers for wall types from original B2E).
    //   u_far_field: The far-field state vector (used if boundary_type_id indicates far-field).
    //   solver: A Riemann solver instance to compute flux for far-field conditions.
    //
    // Returns:
    //   A pair containing:
    //     - The flux vector (4 components) across the boundary face.
    //     - The maximum wave speed associated with this boundary interaction.
    std::pair<std::vector<double>, double> applyBoundaryCondition(
        const std::vector<double>& U_interior,
        const std::vector<double>& normal,
        int boundary_type_id,
        const std::vector<double>& u_far_field, // Provided for far-field conditions
        Solver& solver                          // Riemann solver for far-field
    );

private:
    // Specific boundary condition methods
    // Inviscid wall condition (adapted from wall.h)
    std::pair<std::vector<double>, double> inviscidWall(
        const std::vector<double>& U_interior,
        const std::vector<double>& normal
    );

    // Far-field condition (will use the provided solver)
    std::pair<std::vector<double>, double> farField(
        const std::vector<double>& U_interior,
        const std::vector<double>& u_far_field_state,
        const std::vector<double>& normal,
        Solver& solver
    );

    double gamma_ = 1.4; // Ratio of specific heats
};

#endif // BOUNDARYCONDITIONS_H
