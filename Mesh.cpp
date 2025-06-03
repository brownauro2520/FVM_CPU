#include "Mesh.h"
#include "IOManager.h" // For IOManager class and MeshData struct
#include <iostream>    // For error messages
#include <stdexcept>   // For std::runtime_error

Mesh::Mesh() {
    // Constructor
}

bool Mesh::loadMesh(const std::string& hdf5_filename) {
    IOManager io_manager;
    try {
        IOManager::MeshData mesh_data = io_manager.readMesh(hdf5_filename);

        // Copy data from MeshData struct to Mesh members
        nodes_ = mesh_data.nodes;
        elementsToNodes_ = mesh_data.elementsToNodes;
        cellAreas_ = mesh_data.cellAreas;

        interiorFacesToNodes_ = mesh_data.interiorFacesToNodes;
        interiorFacesToElements_ = mesh_data.interiorFacesToElements;
        interiorFaceNormals_ = mesh_data.interiorFaceNormals;
        interiorFaceLengths_ = mesh_data.interiorFaceLengths;

        boundaryFacesToNodes_ = mesh_data.boundaryFacesToNodes;
        boundaryFacesToElements_ = mesh_data.boundaryFacesToElements;
        boundaryFaceNormals_ = mesh_data.boundaryFaceNormals;
        boundaryFaceLengths_ = mesh_data.boundaryFaceLengths;

        // Optional: Basic validation of loaded data
        if (nodes_.empty() || elementsToNodes_.empty()) {
            std::cerr << "Error: Loaded mesh data appears to be empty or invalid (nodes or elementsToNodes)." << std::endl;
            return false;
        }
        if (elementsToNodes_.size() != cellAreas_.size() && !cellAreas_.empty()) { // cellAreas can be empty if not provided
             std::cerr << "Warning: Number of elements (" << elementsToNodes_.size()
                       << ") does not match number of cell areas (" << cellAreas_.size() << ")." << std::endl;
             // Decide if this is a fatal error or just a warning.
        }
        // Further checks can be added here (e.g., consistency of sizes for faces)
        if (interiorFacesToNodes_.size() != interiorFacesToElements_.size() ||
            interiorFacesToNodes_.size() != interiorFaceNormals_.size() ||
            interiorFacesToNodes_.size() != interiorFaceLengths_.size()) {
            std::cerr << "Warning: Inconsistent sizes for interior face data arrays." << std::endl;
        }
        if (boundaryFacesToNodes_.size() != boundaryFacesToElements_.size() ||
            boundaryFacesToNodes_.size() != boundaryFaceNormals_.size() ||
            boundaryFacesToNodes_.size() != boundaryFaceLengths_.size()) {
            std::cerr << "Warning: Inconsistent sizes for boundary face data arrays." << std::endl;
        }


        std::cout << "Mesh data loaded successfully from (simulated by IOManager) " << hdf5_filename << std::endl;
        std::cout << "Nodes: " << getNumNodes() << ", Elements: " << getNumElements() << std::endl;
        std::cout << "Interior Faces: " << getNumInteriorFaces() << ", Boundary Faces: " << getNumBoundaryFaces() << std::endl;
        if (!cellAreas_.empty()) {
            std::cout << "Cell areas available: Yes (count: " << cellAreas_.size() << ")" << std::endl;
        } else {
            std::cout << "Cell areas available: No" << std::endl;
        }


        return true;
    } catch (const std::runtime_error& e) {
        std::cerr << "Error during mesh loading: " << e.what() << std::endl;
        return false;
    }
}

// Getter implementations
const std::vector<std::vector<double>>& Mesh::getNodes() const { return nodes_; }
const std::vector<std::vector<int>>& Mesh::getElementsToNodes() const { return elementsToNodes_; }
const std::vector<double>& Mesh::getCellAreas() const { return cellAreas_; }
const std::vector<std::vector<int>>& Mesh::getInteriorFacesToNodes() const { return interiorFacesToNodes_; }
const std::vector<std::vector<int>>& Mesh::getInteriorFacesToElements() const { return interiorFacesToElements_; }
const std::vector<std::vector<double>>& Mesh::getInteriorFaceNormals() const { return interiorFaceNormals_; }
const std::vector<double>& Mesh::getInteriorFaceLengths() const { return interiorFaceLengths_; }
const std::vector<std::vector<int>>& Mesh::getBoundaryFacesToNodes() const { return boundaryFacesToNodes_; }
const std::vector<std::vector<int>>& Mesh::getBoundaryFacesToElements() const { return boundaryFacesToElements_; }
const std::vector<std::vector<double>>& Mesh::getBoundaryFaceNormals() const { return boundaryFaceNormals_; }
const std::vector<double>& Mesh::getBoundaryFaceLengths() const { return boundaryFaceLengths_; }

int Mesh::getNumNodes() const { return static_cast<int>(nodes_.size()); }
int Mesh::getNumElements() const { return static_cast<int>(elementsToNodes_.size()); }
// Assuming one entry per face in these vectors
int Mesh::getNumInteriorFaces() const { return static_cast<int>(interiorFacesToNodes_.size()); }
int Mesh::getNumBoundaryFaces() const { return static_cast<int>(boundaryFacesToNodes_.size()); }


// Placeholder geometric calculations (can remain as previously defined or be updated)
double Mesh::calculateCellArea(int elementIndex) const {
    if (elementIndex >= 0 && static_cast<size_t>(elementIndex) < cellAreas_.size()) {
        return cellAreas_[elementIndex]; // Use precomputed if available
    }
    // This part of the function is more of a fallback or for meshes without precomputed areas.
    // Given the current design, cellAreas_ should be populated by IOManager.
    if (elementIndex < 0 || static_cast<size_t>(elementIndex) >= elementsToNodes_.size()) {
        std::cerr << "Error: Invalid element index for area calculation: " << elementIndex << std::endl;
        return 0.0;
    }
    std::cout << "Warning: calculateCellArea called for element " << elementIndex
              << ", but precomputed area was not found or index out of bounds for precomputed areas. Using placeholder." << std::endl;
    return 1.0; // Dummy value
}

std::vector<double> Mesh::calculateFaceNormal(int faceIndex) const {
    // This method is less likely to be used if normals are precomputed and loaded from HDF5.
    // It would need to distinguish between interior and boundary faces if 'faceIndex' is generic
    // and then find the corresponding nodes to calculate the normal.
    // This is complex and error-prone if not done carefully.
    // For now, it's better to rely on normals from IOManager.
    std::cerr << "Error: calculateFaceNormal by generic index is not robustly implemented. "
              << "Rely on precomputed normals (e.g., getInteriorFaceNormals, getBoundaryFaceNormals)." << std::endl;
    return {0.0, 0.0}; // Dummy value
}
