#ifndef FLOWFIELD_H
#define FLOWFIELD_H

#include <vector>
#include <string> // Required for std::string

class FlowField {
public:
    // Constructor: Takes the number of cells and number of variables per cell
    FlowField(int numCells, int numVars);

    // Method to initialize the flow field with uniform values
    void initialize(const std::vector<double>& initialValues);

    // Method to get solution data for a specific cell
    const std::vector<double>& getCellData(int cellIndex) const;

    // Method to get all solution data
    const std::vector<std::vector<double>>& getAllData() const;

    // Method to update solution data for a specific cell
    void updateCellData(int cellIndex, const std::vector<double>& newData);

    // Method to update all solution data
    void updateAllData(const std::vector<std::vector<double>>& newData);

    // Get the number of cells
    int getNumCells() const;

    // Get the number of variables
    int getNumVars() const;


private:
    int numCells_;
    int numVars_;
    std::vector<std::vector<double>> solutionData_; // Stores solution variables for each cell
};

#endif // FLOWFIELD_H
