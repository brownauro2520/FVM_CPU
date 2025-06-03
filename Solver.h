#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <string> // Required by other headers, good to include
#include <utility> // Required for std::pair

// Forward declaration for Mesh and FlowField if needed for method signatures
// class Mesh;
// class FlowField;

class Solver {
public:
    virtual ~Solver() = default; // Important for proper cleanup of derived classes

    // Abstract method for calculating flux
    // UL: Left state vector (rho, rho*u, rho*v, E)
    // UR: Right state vector
    // normal: Normal vector of the face
    // Returns a pair: flux vector and maximum wave speed (for CFL calculation)
    virtual std::pair<std::vector<double>, double> calculateFlux(
        const std::vector<double>& UL,
        const std::vector<double>& UR,
        const std::vector<double>& normal
    ) = 0;

    // Method to get solver type (optional, for debugging or factory patterns)
    virtual std::string getType() const = 0;

protected:
    double gamma_ = 1.4; // Ratio of specific heats, common for both solvers
};

#endif // SOLVER_H
