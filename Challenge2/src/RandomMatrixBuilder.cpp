#include "Matrix.hpp"
#include "MatrixImplementation.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    // Check if the program is being run correctly
    if (argc != 5) {
        std::cout << "Usage : ./matrixBuilder <filename> <rows> <cols> <nnz>" << std::endl;
        return -1;
    }
    // Extract variables inserted through command line
    std::string filename;
    std::size_t num_rows, num_cols, NonZeroElements;
    filename = argv[1];
    num_rows = std::stoul(argv[2]);
    num_cols = std::stoul(argv[3]);
    NonZeroElements = std::stoul(argv[4]);
    // Check that the input data makes sense
    if (num_rows*num_cols < NonZeroElements) {
        std::cout << "Matrix can't store this many elements, run the program with the correct dimensions" << std::endl;
        return -1;
    }
    // Create the empty matrix
    // We use ROW_MAJOR ordering since it's the one that should be used in .mtx files
    algebra::Matrix<double, algebra::ROW_MAJOR> M;
    M.setRows(num_rows);
    M.setCols(num_cols);
    // Start filling the matrix with random values
    std::size_t currRow, currCol;
    double currValue;
    while (M.getNnz() < NonZeroElements) {
        // Generate random indices for value insertion
        currRow = static_cast<std::size_t>(rand() % num_rows);
        currCol = static_cast<std::size_t>(rand() % num_cols);
        // Generate the value to insert
        currValue = static_cast<double>(rand()) / RAND_MAX;
        // Insert the value
        M.set(currRow, currCol, currValue);
    }

    // Export the matrix to .mtx file
    M.SellToMarket(filename);

    return 0;
}