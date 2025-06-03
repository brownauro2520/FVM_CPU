#ifndef HLLESOLVER_H
#define HLLESOLVER_H

#include "Solver.h"
#include <cmath> // For std::sqrt, std::abs, std::max, std::min
#include <algorithm> // For std::max/min

class HLLESolver : public Solver {
public:
    HLLESolver();
    ~HLLESolver() override = default;

    std::pair<std::vector<double>, double> calculateFlux(
        const std::vector<double>& UL,
        const std::vector<double>& UR,
        const std::vector<double>& normal
    ) override;

    std::string getType() const override;
};

#endif // HLLESOLVER_H
