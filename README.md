# 2D Unstructured Finite Volume Method (FVM) CFD Solver (OOP Refactor)

## Objective
This project is an object-oriented C++ refactoring of a 2D unstructured Finite Volume Method (FVM) code for solving compressible Euler equations. The primary goal is to create a modular, extensible, and maintainable framework.

## File Structure

The codebase is organized into several classes, each with a header (`.h`) and an implementation (`.cpp`) file:

*   **`main.cpp`**: The main driver program that sets up and runs the simulation.
*   **`Mesh.h / Mesh.cpp`**: Defines the `Mesh` class, responsible for storing and providing access to mesh data (nodes, element connectivity, face information, areas, normals, etc.).
*   **`FlowField.h / FlowField.cpp`**: Defines the `FlowField` class, which stores and manages the solution variables (density, momentum components, energy) for each cell.
*   **`IOManager.h / IOManager.cpp`**: Defines the `IOManager` class, responsible for mesh input and solution output.
    *   **Note:** Currently, `IOManager::readMesh()` provides hardcoded dummy mesh data. HDF5 reading/writing logic is stubbed and needs to be implemented.
*   **`Solver.h`**: Abstract base class for Riemann solvers.
*   **`RoeSolver.h / RoeSolver.cpp`**: Concrete implementation of the Roe approximate Riemann solver, inheriting from `Solver`.
*   **`HLLESolver.h / HLLESolver.cpp`**: Concrete implementation of the HLLE approximate Riemann solver, inheriting from `Solver`.
*   **`TimeIntegration.h`**: Abstract base class for time-stepping schemes.
*   **`RungeKutta.h / RungeKutta.cpp`**: Concrete implementation of a 4th-order Runge-Kutta time integration scheme, inheriting from `TimeIntegration`.
*   **`BoundaryConditions.h / BoundaryConditions.cpp`**: Defines the `BoundaryConditions` class, responsible for applying boundary conditions (e.g., far-field, inviscid wall) and calculating boundary fluxes.
*   **`extrafunctions.h`**: (Original file) Contains utility functions. Some might be deprecated or integrated into classes.
*   **`readfile.h`**: (Original file) Contains file reading utilities. Likely superseded by `IOManager`.
*   **`wall.h`**: (Original file) Contained inviscid wall logic, now part of `BoundaryConditions`.
*   **`roe.h`, `HLLE.h`**: (Original files) Contained solver logic, now part_of `RoeSolver` and `HLLESolver` respectively.
*   **`residual.h`**: (Original file) Contained residual calculation logic, now integrated into `RungeKutta::calculateResidual`.
*   **`1st_order.cpp`**: (Original main file) This has been replaced by `main.cpp`. It might still be present in the repository if not manually deleted.

## Building the Code

### Dependencies
*   A C++17 compliant compiler (e.g., g++, clang++).
*   **Future:** HDF5 library (when HDF5 I/O is fully implemented in `IOManager`).

### Compilation Example
You can compile the code using a C++ compiler. Ensure all `.cpp` files are included in the compilation command.

```bash
# Example using g++
g++ -std=c++17 -I. -o euler_solver main.cpp Mesh.cpp FlowField.cpp IOManager.cpp RoeSolver.cpp HLLESolver.cpp RungeKutta.cpp BoundaryConditions.cpp
```
*(Note: The `-I.` flag tells the compiler to look for headers in the current directory. Adjust if your headers are elsewhere. The list of .cpp files should include all implemented classes.)*

## Running the Code

Execute the compiled program:

```bash
./euler_solver [CFL_NUMBER] [MAX_ITERATIONS] [SOLVER_CHOICE]
```

### Command-Line Arguments:
*   **`CFL_NUMBER`** (optional, double): The CFL number for time step calculation. Default: `1.2`.
*   **`MAX_ITERATIONS`** (optional, int): The maximum number of iterations to run. Default: `1000`.
*   **`SOLVER_CHOICE`** (optional, int): Selects the Riemann solver.
    *   `0`: Roe solver (default)
    *   `1`: HLLE solver

Example:
```bash
./euler_solver 0.8 5000 0
```
This runs the solver with a CFL of 0.8, for a maximum of 5000 iterations, using the Roe solver.

## Extending the Codebase

The object-oriented design facilitates extensions:

### Adding a New Solver
1.  Create a new class (e.g., `MyNewSolver`) that inherits from the abstract `Solver` class (`Solver.h`).
2.  Implement the constructor and any specific helper methods.
3.  Override the pure virtual function `std::pair<std::vector<double>, double> calculateFlux(...)` to implement your solver's flux calculation logic.
4.  Override `std::string getType() const` to return the name of your solver.
5.  In `main.cpp`, you can then instantiate your new solver: `solver = std::make_unique<MyNewSolver>();`.

### Adding a New Time Integration Scheme
1.  Create a new class (e.g., `MyTimeScheme`) that inherits from `TimeIntegration` (`TimeIntegration.h`).
2.  Implement the constructor.
3.  Override the pure virtual function `double advance(...)` to implement the time-stepping algorithm. This method will typically involve calls to `calculateResidual` (which might become a protected method in `TimeIntegration` or part of a shared utility if more schemes need it directly).
4.  Override `std::string getType() const`.
5.  In `main.cpp`, instantiate your new time integration scheme.

### Adding a New Boundary Condition
1.  Open `BoundaryConditions.h` and `BoundaryConditions.cpp`.
2.  Add a new private method for your boundary condition logic (e.g., `myNewBC(...)`).
3.  Assign a new unique `boundary_type_id`.
4.  Modify the `applyBoundaryCondition` method to dispatch to your new BC method when its `boundary_type_id` is encountered. This ID would typically come from the mesh file.

## Current Limitations
*   **Dummy Mesh Data:** The `IOManager` currently provides a hardcoded, simple dummy mesh. Functionality to read mesh data from external files (e.g., HDF5) is not yet implemented.
*   **Stubbed Solution Output:** Solution writing to HDF5 is also stubbed in `IOManager`.
*   **2D Only:** The framework is designed for 2D, though 3D extension was a design consideration.
*   **No Parallelization:** The code is serial. Parallelization (MPI, OpenMP) is a future extension goal.

```
