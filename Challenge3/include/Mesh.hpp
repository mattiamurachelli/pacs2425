#ifndef MESH_HPP
#define MESH_HPP

#include <utility>
#include <vector>

// Class that represents a square mesh of side L=1
class Mesh {
    public :
        using coords_type = std::pair<double,double>;
        using coordinates_matrix = std::vector<std::vector<coords_type>>;

        // Constructor
        Mesh(std::size_t size);
        // Getter for mesh coordinates
        coords_type operator() (std::size_t i, std::size_t j) const;
        // Getter for mesh size
        double get_h() const;

    private :
        // Data
        std::size_t size;                                       // Number of rows and columns of the matrix
        double h;                                               // Mesh size
        coordinates_matrix physical_points;                     // Coordinate points in the mesh
};

#endif // MESH_HPP