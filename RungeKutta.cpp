#include "RungeKutta.h"
#include "FlowField.h"
#include "Mesh.h"
#include "Solver.h"
#include "BoundaryConditions.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <algorithm> // For std::floor

// Static helper functions (can be defined in .cpp file)
static std::vector<double> flowFieldToVector(const FlowField& ff) {
    int nCells = ff.getNumCells();
    int nVars = ff.getNumVars();
    std::vector<double> vec(nCells * nVars);
    const auto& data = ff.getAllData();
    for (int i = 0; i < nCells; ++i) {
        for (int j = 0; j < nVars; ++j) {
            vec[i * nVars + j] = data[i][j];
        }
    }
    return vec;
}

static void vectorToFlowField(const std::vector<double>& vec, FlowField& ff) {
    int nCells = ff.getNumCells();
    int nVars = ff.getNumVars();
    if (vec.size() != static_cast<size_t>(nCells * nVars)) { // Explicit cast for comparison
        throw std::runtime_error("Vector size does not match FlowField dimensions in vectorToFlowField. Vector size: " + std::to_string(vec.size()) + ", FF dims: " + std::to_string(nCells) + "x" + std::to_string(nVars));
    }
    std::vector<std::vector<double>> data(nCells, std::vector<double>(nVars));
    for (int i = 0; i < nCells; ++i) {
        for (int j = 0; j < nVars; ++j) {
            data[i][j] = vec[static_cast<size_t>(i) * nVars + j];
        }
    }
    ff.updateAllData(data);
}


RungeKutta::RungeKutta(int order) : order_(order) {
    if (order_ != 4) {
        throw std::invalid_argument("RungeKutta currently only supports 4th order.");
    }
}

std::string RungeKutta::getType() const {
    return "RungeKutta" + std::to_string(order_);
}

void RungeKutta::calculateResidual(
    std::vector<double>& R_out,
    std::vector<double>& dT_out,
    const FlowField& U_state_ff,
    const Mesh& mesh_obj,
    Solver& solver,
    BoundaryConditions& bc,
    const std::vector<double>& u_inf_state,
    double CFL_val) {

    int N_cells = U_state_ff.getNumCells();
    int N_vars = U_state_ff.getNumVars();

    if (N_vars != 4) {
        throw std::runtime_error("Residual calculation expects 4 variables per cell.");
    }
    if (u_inf_state.size() != 4 && !u_inf_state.empty()) { // u_inf_state can be empty if no far-field BCs
        // This check might be too strict if u_inf_state is passed but not used by any active BC.
        // However, if a far-field BC is hit, BoundaryConditions class will throw error if it's wrong size.
        // std::cout << "Warning: u_inf_state size is not 4 in calculateResidual. Size: " << u_inf_state.size() << std::endl;
    }


    R_out.assign(N_cells * N_vars, 0.0);
    dT_out.assign(N_cells, 0.0);

    const auto& I2E = mesh_obj.getInteriorFacesToElements();
    const auto& interior_normals = mesh_obj.getInteriorFaceNormals();
    const auto& interior_lengths = mesh_obj.getInteriorFaceLengths();

    for (size_t i = 0; i < static_cast<size_t>(mesh_obj.getNumInteriorFaces()); ++i) {
        if (I2E.size() <= i || I2E[i].size() < 2) {
            throw std::runtime_error("Interior face to elements data (I2E) is incomplete or missing for face " + std::to_string(i));
        }
        int ul_idx = I2E[i][0];
        int ur_idx = I2E[i][1];

        if (ul_idx < 0 || ul_idx >= N_cells || ur_idx < 0 || ur_idx >= N_cells) {
            throw std::out_of_range("Element index out of range for interior face " + std::to_string(i) +
                                    ". ul_idx: " + std::to_string(ul_idx) + ", ur_idx: " + std::to_string(ur_idx) +
                                    ", N_cells: " + std::to_string(N_cells));
        }

        const std::vector<double>& UL_vars = U_state_ff.getCellData(ul_idx);
        const std::vector<double>& UR_vars = U_state_ff.getCellData(ur_idx);

        if (interior_normals.size() <= i || interior_normals[i].size() < 2) {
             throw std::runtime_error("Interior face normal data is incomplete or missing for face " + std::to_string(i));
        }
        const std::vector<double>& normal_vec = interior_normals[i];

        if (interior_lengths.size() <= i) {
            throw std::runtime_error("Interior face length data is missing for face " + std::to_string(i));
        }
        double edge_length = interior_lengths[i];
        if (edge_length <= 0) {
            throw std::runtime_error("Non-positive edge length for interior face " + std::to_string(i));
        }


        std::pair<std::vector<double>, double> flux_result = solver.calculateFlux(UL_vars, UR_vars, normal_vec);
        const auto& F_flux = flux_result.first;
        double s_max = flux_result.second;

        for (int j = 0; j < N_vars; ++j) {
            R_out[static_cast<size_t>(ul_idx) * N_vars + j] += F_flux[j] * edge_length;
            R_out[static_cast<size_t>(ur_idx) * N_vars + j] -= F_flux[j] * edge_length;
        }
        dT_out[ul_idx] += std::abs(s_max) * edge_length;
        dT_out[ur_idx] += std::abs(s_max) * edge_length;
    }

    const auto& B2E = mesh_obj.getBoundaryFacesToElements();
    const auto& boundary_normals = mesh_obj.getBoundaryFaceNormals();
    const auto& boundary_lengths = mesh_obj.getBoundaryFaceLengths();

    for (size_t i = 0; i < static_cast<size_t>(mesh_obj.getNumBoundaryFaces()); ++i) {
        if (B2E.size() <= i || B2E[i].size() < 2) {
            throw std::runtime_error("Boundary face to elements data (B2E) is incomplete or missing for face " + std::to_string(i));
        }
        int ul_idx = B2E[i][0];
        int boundary_type_id = B2E[i][1];

        if (ul_idx < 0 || ul_idx >= N_cells) {
             throw std::out_of_range("Element index out of range for boundary face " + std::to_string(i) +
                                     ". ul_idx: " + std::to_string(ul_idx) + ", N_cells: " + std::to_string(N_cells));
        }

        const std::vector<double>& UL_vars = U_state_ff.getCellData(ul_idx);
        if (boundary_normals.size() <= i || boundary_normals[i].size() < 2) {
            throw std::runtime_error("Boundary face normal data is incomplete or missing for face " + std::to_string(i));
        }
        const std::vector<double>& normal_vec = boundary_normals[i];
        if (boundary_lengths.size() <= i) {
            throw std::runtime_error("Boundary face length data is missing for face " + std::to_string(i));
        }
        double edge_length = boundary_lengths[i];
        if (edge_length <= 0) {
            throw std::runtime_error("Non-positive edge length for boundary face " + std::to_string(i));
        }

        std::pair<std::vector<double>, double> flux_result = bc.applyBoundaryCondition(
            UL_vars, normal_vec, boundary_type_id, u_inf_state, solver
        );
        const auto& F_flux = flux_result.first;
        double s_max = flux_result.second;

        for (int j = 0; j < N_vars; ++j) {
            R_out[static_cast<size_t>(ul_idx) * N_vars + j] += F_flux[j] * edge_length;
        }
        dT_out[ul_idx] += std::abs(s_max) * edge_length;
    }

    const auto& cell_areas = mesh_obj.getCellAreas();
    if (cell_areas.size() != static_cast<size_t>(N_cells)) {
        throw std::runtime_error("Mismatch between number of cells (" + std::to_string(N_cells) +
                                 ") and cell areas data size (" + std::to_string(cell_areas.size()) + ").");
    }

    for (int i = 0; i < N_cells; ++i) {
        if (dT_out[i] < 1e-12) { // Sum of s_max * L is very small or zero
            // This cell might be isolated or have zero wave speeds on all faces.
            // Set a very large dt to effectively freeze this cell's state for this step.
            // Or throw error if this is unexpected.
            // std::cerr << "Warning: Sum of s_max*L is near zero for cell " << i << ". Setting dt to large value." << std::endl;
            dT_out[i] = 1e12;
        } else {
            double area = cell_areas[i];
            if (area <= 1e-12) {
                throw std::runtime_error("Non-positive or zero cell area encountered for cell " + std::to_string(i) + ": " + std::to_string(area));
            }
            // The original code used 2.0 * CFL. Standard CFL definition is often CFL * Area / sum(s_max*L)
            // For this exercise, using 2.0 * CFL to match the spirit of the original code.
            dT_out[i] = (2.0 * CFL_val * area) / dT_out[i];
        }
    }
}

double RungeKutta::advance(
    FlowField& U_ff,
    const Mesh& mesh,
    Solver& solver,
    BoundaryConditions& bc,
    const std::vector<double>& u_inf_state,
    double CFL,
    double& totalResidualStore
) {
    int N_cells = U_ff.getNumCells();
    int N_vars = U_ff.getNumVars();

    if (N_vars != 4) {
        throw std::runtime_error("RungeKutta advance expects 4 variables per cell in FlowField.");
    }

    std::vector<double> U_vec_initial = flowFieldToVector(U_ff);
    std::vector<double> k1(static_cast<size_t>(N_cells) * N_vars), k2(static_cast<size_t>(N_cells) * N_vars), k3(static_cast<size_t>(N_cells) * N_vars), k4(static_cast<size_t>(N_cells) * N_vars);
    std::vector<double> u_temp_vec(static_cast<size_t>(N_cells) * N_vars);
    std::vector<double> R_stage(static_cast<size_t>(N_cells) * N_vars);
    std::vector<double> dt_for_stage(N_cells);

    FlowField U_temp_ff(N_cells, N_vars);

    const auto& cell_areas = mesh.getCellAreas();
    if (cell_areas.size() != static_cast<size_t>(N_cells)) {
        throw std::runtime_error("Mismatch between number of cells and cell areas data in advance(). Loaded areas: " + std::to_string(cell_areas.size()) + ", N_cells: " + std::to_string(N_cells) );
    }

    // Stage 1: Calculate k1
    vectorToFlowField(U_vec_initial, U_temp_ff);
    calculateResidual(R_stage, dt_for_stage, U_temp_ff, mesh, solver, bc, u_inf_state, CFL);
    for (int i = 0; i < N_cells * N_vars; ++i) {
        int cell_idx = static_cast<int>(std::floor(static_cast<double>(i) / N_vars));
        double area = cell_areas[cell_idx];
        if (area <= 1e-12) throw std::runtime_error("Non-positive cell area for cell " + std::to_string(cell_idx) + " in k1 calc.");
        k1[i] = R_stage[i] / (-area); // k = dU/dt = -Res/Area
    }

    // Stage 2: Calculate k2
    for (int i = 0; i < N_cells * N_vars; ++i) {
        int cell_idx = static_cast<int>(std::floor(static_cast<double>(i) / N_vars));
        u_temp_vec[i] = U_vec_initial[i] + 0.5 * k1[i] * dt_for_stage[cell_idx];
    }
    vectorToFlowField(u_temp_vec, U_temp_ff);
    // Note: dt_for_stage is typically calculated once based on the initial state U_vec_initial for a standard RK4 step.
    // Re-calling calculateResidual will re-calculate dt_for_stage based on U_temp_ff.
    // For simplicity here, we allow dt_for_stage to be recomputed, though for classic RK, it's fixed.
    // If a fixed dt per cell for the entire RK step is desired, the first dt_for_stage should be stored and reused.
    calculateResidual(R_stage, dt_for_stage, U_temp_ff, mesh, solver, bc, u_inf_state, CFL);
    for (int i = 0; i < N_cells * N_vars; ++i) {
        int cell_idx = static_cast<int>(std::floor(static_cast<double>(i) / N_vars));
        double area = cell_areas[cell_idx];
        if (area <= 1e-12) throw std::runtime_error("Non-positive cell area for cell " + std::to_string(cell_idx) + " in k2 calc.");
        k2[i] = R_stage[i] / (-area);
    }

    // Stage 3: Calculate k3
    for (int i = 0; i < N_cells * N_vars; ++i) {
        int cell_idx = static_cast<int>(std::floor(static_cast<double>(i) / N_vars));
        u_temp_vec[i] = U_vec_initial[i] + 0.5 * k2[i] * dt_for_stage[cell_idx];
    }
    vectorToFlowField(u_temp_vec, U_temp_ff);
    calculateResidual(R_stage, dt_for_stage, U_temp_ff, mesh, solver, bc, u_inf_state, CFL);
    for (int i = 0; i < N_cells * N_vars; ++i) {
        int cell_idx = static_cast<int>(std::floor(static_cast<double>(i) / N_vars));
        double area = cell_areas[cell_idx];
        if (area <= 1e-12) throw std::runtime_error("Non-positive cell area for cell " + std::to_string(cell_idx) + " in k3 calc.");
        k3[i] = R_stage[i] / (-area);
    }

    // Stage 4: Calculate k4
    for (int i = 0; i < N_cells * N_vars; ++i) {
        int cell_idx = static_cast<int>(std::floor(static_cast<double>(i) / N_vars));
        u_temp_vec[i] = U_vec_initial[i] + k3[i] * dt_for_stage[cell_idx];
    }
    vectorToFlowField(u_temp_vec, U_temp_ff);
    calculateResidual(R_stage, dt_for_stage, U_temp_ff, mesh, solver, bc, u_inf_state, CFL);
    for (int i = 0; i < N_cells * N_vars; ++i) {
        int cell_idx = static_cast<int>(std::floor(static_cast<double>(i) / N_vars));
        double area = cell_areas[cell_idx];
        if (area <= 1e-12) throw std::runtime_error("Non-positive cell area for cell " + std::to_string(cell_idx) + " in k4 calc.");
        k4[i] = R_stage[i] / (-area);
    }

    std::vector<double> U_vec_final = U_vec_initial;
    for (int i = 0; i < N_cells * N_vars; ++i) {
        int cell_idx = static_cast<int>(std::floor(static_cast<double>(i) / N_vars));
        U_vec_final[i] += (1.0 / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * dt_for_stage[cell_idx];
    }
    vectorToFlowField(U_vec_final, U_ff);

    totalResidualStore = 0.0;
    // The R_stage from the last calculateResidual call (for k4) is typically used
    // to assess convergence or as the residual for the time step.
    for (size_t i = 0; i < static_cast<size_t>(N_cells * N_vars); ++i) {
        totalResidualStore += std::abs(R_stage[i]);
    }

    return totalResidualStore;
}
