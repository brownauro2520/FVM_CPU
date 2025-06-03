#ifndef IOMANAGER_H
#define IOMANAGER_H

#include <string>
#include <vector>
#include <map> // For returning structured mesh data

// Forward declarations
class FlowField;
class Mesh; // Mesh class will be populated by IOManager

// Define a structure to hold all mesh data components read from HDF5
// This structure will be used to pass data to the Mesh class.
struct MeshData {
    // Basic geometry
    std::vector<std::vector<double>> nodes;       // Nx2 or Nx3 for 2D/3D
    std::vector<std::vector<int>> elementsToNodes; // Element connectivity (e.g., triangles, quads)
    std::vector<double> cellAreas;

    // Interior faces/edges
    std::vector<std::vector<int>> interiorFacesToNodes; // Connectivity of interior faces
    std::vector<std::vector<int>> interiorFacesToElements; // Left/Right elements for each interior face
    std::vector<std::vector<double>> interiorFaceNormals;
    std::vector<double> interiorFaceLengths;

    // Boundary faces/edges
    std::vector<std::vector<int>> boundaryFacesToNodes; // Connectivity of boundary faces
    std::vector<std::vector<int>> boundaryFacesToElements; // Parent element and BC type ID
    std::vector<std::vector<double>> boundaryFaceNormals;
    std::vector<double> boundaryFaceLengths;

    // Add any other mesh components if necessary
    // For example, node-to-element connectivity, etc.
};


class IOManager {
public:
    IOManager();
    ~IOManager() = default;

    // Reads mesh data from an HDF5 file and returns it in a MeshData struct.
    // The Mesh object will then use this struct to populate itself.
    // Throws std::runtime_error on failure.
    MeshData readMesh(const std::string& hdf5_filename);

    // Writes solution data from FlowField to an HDF5 file.
    // timestep_group_name could be e.g., "Timestep_00100"
    // Returns true on success, false on failure.
    bool writeSolution(
        const std::string& hdf5_filename,
        const FlowField& flow_field,
        const std::string& timestep_group_name // e.g., "Solution/Timestep_00100"
    );

    // Helper: Check if an HDF5 file exists and is valid
    bool checkFile(const std::string& hdf5_filename);

private:
    // Private helper methods for HDF5 operations.
    // These would use the HDF5 C API.
    // Example:
    // std::vector<std::vector<double>> readDoubleDataset(const std::string& filename, const std::string& dataset_path);
    // std::vector<double> readDouble1DDataset(const std::string& filename, const std::string& dataset_path);
    // std::vector<std::vector<int>> readIntDataset(const std::string& filename, const std::string& dataset_path);
    // bool writeDoubleDataset(const std::string& filename, const std::string& dataset_path,
    //                         const std::vector<std::vector<double>>& data,
    //                         const std::vector<hsize_t>& dims);

    // For this subtask, the actual HDF5 calls will be stubbed out or simplified
    // to focus on the class structure and method signatures.
    // A full implementation would require linking against an HDF5 library.
};

#endif // IOMANAGER_H
