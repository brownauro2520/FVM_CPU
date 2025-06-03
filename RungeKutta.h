#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "TimeIntegration.h"
#include <vector>

// Forward declarations
class FlowField;
class Mesh;
class Solver;
class BoundaryConditions;

class RungeKutta : public TimeIntegration {
public:
    RungeKutta(int order = 4);
    ~RungeKutta() override = default;

    double advance(
        FlowField& U,
        const Mesh& mesh,
        Solver& solver,
        BoundaryConditions& bc,
        const std::vector<double>& u_inf_state, // MODIFIED: Added u_inf_state
        double CFL,
        double& totalResidualStore
    ) override;

    std::string getType() const override;

private:
    int order_;

    void calculateResidual(
        std::vector<double>& R_out,
        std::vector<double>& dT_out,
        const FlowField& U_state_ff,
        const Mesh& mesh_obj,
        Solver& solver,
        BoundaryConditions& bc,
        const std::vector<double>& u_inf_state,
        double CFL_val
    );
};

#endif // RUNGEKUTTA_H
