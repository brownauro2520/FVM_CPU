#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include <memory>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "Mesh.h"
#include "FlowField.h"
#include "Solver.h"
#include "RoeSolver.h"
#include "HLLESolver.h"
#include "TimeIntegration.h"
#include "RungeKutta.h"
#include "BoundaryConditions.h"
#include "IOManager.h"


int main(int argc, char* argv[]) {
    std::cout << "Starting 2D Compressible Euler Solver (OOP Refactor)" << std::endl;

    double mach_number = 0.25;
    double angle_of_attack_deg = 0.0;
    double cfl_number = 1.2;
    double gamma = 1.4;
    int max_iterations = 1000; // Reduced for quick testing with dummy data
    double convergence_tolerance = 1e-5;
    std::string mesh_filename = "mesh_data.h5"; // Will use dummy data from IOManager
    std::string output_filename = "solution_output.h5";
    int output_frequency = 100;
    int solver_choice = 0; // 0 for Roe, 1 for HLLE

    if (argc > 1) {
        try {
            cfl_number = std::stod(argv[1]);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not parse CFL from command line argument '" << argv[1] << "'. Using default: " << cfl_number << std::endl;
        }
    }
    if (argc > 2) {
        try {
            max_iterations = std::stoi(argv[2]);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not parse Max Iterations from command line argument '" << argv[2] << "'. Using default: " << max_iterations << std::endl;
        }
    }
     if (argc > 3) {
        try {
            solver_choice = std::stoi(argv[3]);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not parse Solver Choice (0=Roe, 1=HLLE) from command line argument '" << argv[3] << "'. Using default: Roe" << std::endl;
        }
    }

    Mesh mesh;
    IOManager io_manager;
    BoundaryConditions boundary_conditions;

    std::unique_ptr<Solver> solver;
    if (solver_choice == 1) {
        solver = std::make_unique<HLLESolver>();
    } else {
        solver = std::make_unique<RoeSolver>();
    }
    std::cout << "Using Solver: " << solver->getType() << std::endl;

    std::unique_ptr<TimeIntegration> time_integrator = std::make_unique<RungeKutta>(4);
    std::cout << "Using Time Integrator: " << time_integrator->getType() << std::endl;

    int final_iter_count = 0;

    try {
        std::cout << "Loading mesh from: " << mesh_filename << "..." << std::endl;
        // IOManager::readMesh() currently returns dummy data, so mesh_filename is a placeholder.
        if (!mesh.loadMesh(mesh_filename)) {
            std::cerr << "Error: Mesh loading failed." << std::endl;
            return 1;
        }
        std::cout << "Mesh loaded successfully." << std::endl;
        std::cout << "  Nodes: " << mesh.getNumNodes() << ", Elements: " << mesh.getNumElements() << std::endl;
        std::cout << "  Interior Faces: " << mesh.getNumInteriorFaces() << ", Boundary Faces: " << mesh.getNumBoundaryFaces() << std::endl;

        if (mesh.getNumElements() == 0) {
            std::cerr << "Error: Mesh contains no elements. Cannot proceed." << std::endl;
            return 1;
        }

        int num_cells = mesh.getNumElements();
        int num_vars = 4;
        FlowField flow_field(num_cells, num_vars);

        double angle_of_attack_rad = angle_of_attack_deg * M_PI / 180.0;
        double gmi = gamma - 1.0;
        if (gmi <= 0 || gamma <=0) { // Ensure gamma is valid before using
            throw std::runtime_error("Gamma value must be greater than 1.");
        }
        std::vector<double> u_inf_state = {
            1.0,
            mach_number * std::cos(angle_of_attack_rad),
            mach_number * std::sin(angle_of_attack_rad),
            (1.0 / (gmi * gamma)) + (0.5 * mach_number * mach_number)
        };

        std::cout << "Initializing flow field with far-field state:" << std::endl;
        std::cout << "  u_inf = [" << u_inf_state[0] << ", " << u_inf_state[1] << ", "
                  << u_inf_state[2] << ", " << u_inf_state[3] << "]" << std::endl;
        flow_field.initialize(u_inf_state);

        std::cout << "\nStarting simulation..." << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << std::fixed << std::setprecision(8);

        auto start_time = std::chrono::high_resolution_clock::now();
        double current_residual_tracker = 0.0; // Use a separate variable for the output parameter of advance

        for (int iter = 0; iter < max_iterations; ++iter) {
            final_iter_count = iter + 1;
            // Pass current_residual_tracker as the output parameter
            time_integrator->advance(
                flow_field,
                mesh,
                *solver,
                boundary_conditions,
                u_inf_state,
                cfl_number,
                current_residual_tracker
            );

            std::cout << "Iter: " << std::setw(5) << iter
                      << " | Residual: " << std::scientific << current_residual_tracker << std::defaultfloat << std::endl;

            // Output solution at specified frequency
            if ((iter + 1) % output_frequency == 0) {
                std::string sol_group_name = "Solution/Timestep_" + std::to_string(iter + 1);
                //std::cout << "Writing solution at iteration " << iter + 1 << " to group " << sol_group_name << std::endl;
                if (!io_manager.writeSolution(output_filename, flow_field, sol_group_name)) { // IOManager::writeSolution is currently stubbed
                    std::cerr << "Warning: Failed to write solution at iteration " << iter + 1 << std::endl;
                }
            }

            // Convergence Check
            if (current_residual_tracker < convergence_tolerance) {
                std::cout << "\nConvergence reached at iteration " << iter << "!" << std::endl;
                 // Write final converged solution if it wasn't on an output frequency
                 if (((iter + 1) % output_frequency != 0) ) {
                    io_manager.writeSolution(output_filename, flow_field, "Solution/FinalConverged_" + std::to_string(iter+1));
                }
                break;
            }

            // Max iterations reached
             if (iter == max_iterations - 1) {
                std::cout << "\nMaximum iterations reached." << std::endl;
                // Write final solution if it wasn't on an output frequency
                if (((iter + 1) % output_frequency != 0) ) {
                     io_manager.writeSolution(output_filename, flow_field, "Solution/FinalMaxIter_" + std::to_string(iter+1));
                }
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Simulation finished." << std::endl;
        std::cout << "Total iterations performed: " << final_iter_count << std::endl;
        std::cout << "Final Residual: " << std::scientific << current_residual_tracker << std::defaultfloat << std::endl;
        std::cout << "Total simulation time: " << std::fixed << std::setprecision(3) << elapsed_time.count() << " seconds." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\n!!! An exception occurred: " << e.what() << " !!!" << std::endl;
        return 1;
    }

    std::cout << "\nSolver execution completed." << std::endl;
    return 0;
}
