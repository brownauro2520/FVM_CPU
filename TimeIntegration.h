#ifndef TIMEINTEGRATION_H
#define TIMEINTEGRATION_H

#include <vector>
#include <string> // For getType

// Forward declarations
class FlowField;
class Mesh;
class Solver; // For residual calculation
class BoundaryConditions; // For residual calculation

class TimeIntegration {
public:
    virtual ~TimeIntegration() = default;

    // Method to advance the solution by one time step
    // U: FlowField object (contains solution data, will be updated)
    // mesh: Mesh object (provides geometry info)
    // solver: Solver object (to calculate fluxes)
    // bc: BoundaryConditions object
    // CFL: CFL number
    // Returns the total residual after the time step
    virtual double advance(
        FlowField& U,
        const Mesh& mesh,
        Solver& solver,
        BoundaryConditions& bc,
        double CFL,
        double& totalResidual // Output parameter for total residual
    ) = 0;

    virtual std::string getType() const = 0;

protected:
    double gamma_ = 1.4; // Might be needed by residual calculation logic
};

#endif // TIMEINTEGRATION_H
