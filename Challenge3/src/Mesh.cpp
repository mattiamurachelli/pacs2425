#include "Mesh.hpp"
#include <iostream>

// Constructor -> generates the mesh
Mesh::Mesh(std::size_t size_) {
    // Assign number of rows and columns
    size = size_;

    // Compute the mesh size
    h = 1.0 / (size-1);

    // Compute the physical points coordinates
    physical_points.resize(size, std::vector<coords_type>(size));
    for(std::size_t i=0; i<size; ++i) {
        for(std::size_t j=0; j<size; ++j) {
            physical_points[i][j].first = i*h;
            physical_points[i][j].second = j*h;
        }
    }
}
// Getter for mesh coordinates
Mesh::coords_type Mesh::operator() (std::size_t i, std::size_t j) const {
    return physical_points[i][j];
}
// Getter for mesh size
double Mesh::get_h() const {
    return h;
}