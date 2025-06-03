#include "FlowField.h"
#include <stdexcept> // Required for std::out_of_range

// Constructor
FlowField::FlowField(int numCells, int numVars) : numCells_(numCells), numVars_(numVars) {
    if (numCells <= 0 || numVars <= 0) {
        throw std::invalid_argument("Number of cells and variables must be positive.");
    }
    solutionData_.resize(numCells, std::vector<double>(numVars, 0.0));
}

// Method to initialize the flow field with uniform values
void FlowField::initialize(const std::vector<double>& initialValues) {
    if (initialValues.size() != numVars_) {
        throw std::invalid_argument("Size of initialValues must match the number of variables.");
    }
    for (int i = 0; i < numCells_; ++i) {
        solutionData_[i] = initialValues;
    }
}

// Method to get solution data for a specific cell
const std::vector<double>& FlowField::getCellData(int cellIndex) const {
    if (cellIndex < 0 || cellIndex >= numCells_) {
        throw std::out_of_range("Cell index out of range.");
    }
    return solutionData_[cellIndex];
}

// Method to get all solution data
const std::vector<std::vector<double>>& FlowField::getAllData() const {
    return solutionData_;
}

// Method to update solution data for a specific cell
void FlowField::updateCellData(int cellIndex, const std::vector<double>& newData) {
    if (cellIndex < 0 || cellIndex >= numCells_) {
        throw std::out_of_range("Cell index out of range.");
    }
    if (newData.size() != numVars_) {
        throw std::invalid_argument("Size of newData must match the number of variables.");
    }
    solutionData_[cellIndex] = newData;
}

// Method to update all solution data
void FlowField::updateAllData(const std::vector<std::vector<double>>& newData) {
    if (newData.size() != numCells_) {
        throw std::invalid_argument("Number of rows in newData must match the number of cells.");
    }
    for(int i = 0; i < numCells_; ++i) {
        if (newData[i].size() != numVars_) {
            throw std::invalid_argument("Number of columns in newData must match the number of variables for cell " + std::to_string(i));
        }
    }
    solutionData_ = newData;
}

// Get the number of cells
int FlowField::getNumCells() const {
    return numCells_;
}

// Get the number of variables
int FlowField::getNumVars() const {
    return numVars_;
}
