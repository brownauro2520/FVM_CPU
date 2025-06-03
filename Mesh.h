#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include "IOManager.h" // For MeshData

class Mesh {
public:
    Mesh();
    bool loadMesh(const std::string& hdf5_filename);

    // Getters
    const std::vector<std::vector<double>>& getNodes() const;
    const std::vector<std::vector<int>>& getElementsToNodes() const;
    const std::vector<double>& getCellAreas() const;
    const std::vector<std::vector<int>>& getInteriorFacesToNodes() const;
    const std::vector<std::vector<int>>& getInteriorFacesToElements() const;
    const std::vector<std::vector<double>>& getInteriorFaceNormals() const;
    const std::vector<double>& getInteriorFaceLengths() const;
    const std::vector<std::vector<int>>& getBoundaryFacesToNodes() const;
    const std::vector<std::vector<int>>& getBoundaryFacesToElements() const;
    const std::vector<std::vector<double>>& getBoundaryFaceNormals() const;
    const std::vector<double>& getBoundaryFaceLengths() const;

    int getNumNodes() const;
    int getNumElements() const;
    int getNumInteriorFaces() const;
    int getNumBoundaryFaces() const;

    // Placeholder geometric calculations (can be refined or removed if data is always precomputed)
    double calculateCellArea(int elementIndex) const; // May use precomputed cellAreas_
    std::vector<double> calculateFaceNormal(int faceIndex) const; // May use precomputed normals

private:
    std::vector<std::vector<double>> nodes_;
    std::vector<std::vector<int>> elementsToNodes_;
    std::vector<double> cellAreas_;

    std::vector<std::vector<int>> interiorFacesToNodes_;
    std::vector<std::vector<int>> interiorFacesToElements_; // Stores {left_elem_idx, right_elem_idx}
    std::vector<std::vector<double>> interiorFaceNormals_;
    std::vector<double> interiorFaceLengths_;

    std::vector<std::vector<int>> boundaryFacesToNodes_;
    std::vector<std::vector<int>> boundaryFacesToElements_; // Stores {elem_idx, boundary_type_id}
    std::vector<std::vector<double>> boundaryFaceNormals_;
    std::vector<double> boundaryFaceLengths_;
};
#endif // MESH_H
