#ifndef ROESOLVER_H
#define ROESOLVER_H

#include "Solver.h"
#include <cmath> // For std::sqrt, std::abs, etc.
#include <algorithm> // For std::max_element, std::min etc.

class RoeSolver : public Solver {
public:
    RoeSolver();
    ~RoeSolver() override = default;

    std::pair<std::vector<double>, double> calculateFlux(
        const std::vector<double>& UL,
        const std::vector<double>& UR,
        const std::vector<double>& normal
    ) override;

    std::string getType() const override;

private:
    // Helper functions from the original roe.h, adapted as private methods
    std::vector<double> absolute_value(const std::vector<double>& v);
    std::vector<double> subtract(const std::vector<double>& v1, const std::vector<double>& v2);
};

#endif // ROESOLVER_H
