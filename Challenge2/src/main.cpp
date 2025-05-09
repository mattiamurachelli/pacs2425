#include "chrono.hpp"
#include "Matrix.hpp"
#include "MatrixImplementation.hpp"

#include <iostream>

std::vector<double> randomVector(std::size_t size) {
    std::vector<double> v(size);
    for (std::size_t i = 0; i < size; ++i) {
        v[i] = static_cast<double>(rand()) / RAND_MAX;
    }
    return v;
}

void printVector(std::vector<double> const &v) {
    for (std::size_t i = 0; i < v.size(); ++i) {
        std::cout << v[i] << "\t";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]){
    // FILE FOR SHOWING THE DIFFERENCE IN PERFORMANCE OF THE RIGHT VECTOR MULTIPLICATION
    // AND NORM FUNCTIONS BETWEEN COMPRESSED AND UNCOMPRESSED MATRICES

    // Get from the user filename and option to print verbose outputs
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage : ./main <filename> Y(/N)" << std::endl;
        return -1;
    }

    // Assign filename and check if file exists
    std::string filename = "data/" + std::string(argv[1]);
    std::ifstream file_check(filename);
    if (!file_check.good()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return -1;
    }

    // Assign output options (verbose or non verbose)
    bool printOutput = false;
    if (argc == 3 && std::string(argv[2]) == "Y"){
        printOutput = true;
    }

    // Start creating a timer object
    Timings::Chrono timer;

    // Row Major Matrix with double coefficients
    std::cout << "=================================================================" << std::endl;
    std::cout << "Row Major Matrix with double coefficients" << std::endl;
    algebra::Matrix<double, algebra::ROW_MAJOR> M(filename);
    if (printOutput) {M.print();}

    // We create a random vector for the multiplication
    auto v = randomVector(M.getCols());
    if (printOutput) {
        std::cout << "Vector for multiplication : " << std::endl;
        printVector(v);
    }

    // And define the double to store norms
    double norm;
    
    // We begin with uncompressed state
    // Perform the multiplication
    std::cout << "Right vector multiplication in uncompressed state" << std::endl;
    timer.start();
    auto result_row_u = M * v;
    timer.stop();
    if (printOutput) {
        std::cout << "Result : " << std::endl;
        printVector(result_row_u);
    }
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Perform norm computations
    // One norm
    timer.start();
    norm = M.norm<algebra::ONE>();
    timer.stop();
    std::cout << "ONE norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Infty norm
    timer.start();
    norm = M.norm<algebra::INFTY>();
    timer.stop();
    std::cout << "INFTY norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Frobenius norm
    timer.start();
    norm = M.norm<algebra::FROBENIUS>();
    timer.stop();
    std::cout << "FROBENIUS norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    std::cout << "=================================================================" << std::endl;

    // Compress the matrix
    M.compress();

    // Perform the multiplication
    std::cout << "Right vector multiplication in compressed state" << std::endl;
    timer.start();
    auto result_row_c = M * v;
    timer.stop();
    if (printOutput) {
        std::cout << "Result : " << std::endl;
        printVector(result_row_c);
    }
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Perform norm computations
    // One norm
    timer.start();
    norm = M.norm<algebra::ONE>();
    timer.stop();
    std::cout << "ONE norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Infty norm
    timer.start();
    norm = M.norm<algebra::INFTY>();
    timer.stop();
    std::cout << "INFTY norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Frobenius norm
    timer.start();
    norm = M.norm<algebra::FROBENIUS>();
    timer.stop();
    std::cout << "FROBENIUS norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    std::cout << "=================================================================" << std::endl;

    /********************************************************************************************/

    // Column Major Matrix with double coefficients
    std::cout << "=================================================================" << std::endl;
    std::cout << "Column Major Matrix with double coefficients" << std::endl;
    algebra::Matrix<double, algebra::COL_MAJOR> N(filename);
    if (printOutput) {N.print();}

    // We begin with uncompressed state
    // Perform the multiplication
    std::cout << "Right vector multiplication in uncompressed state" << std::endl;
    timer.start();
    auto result_col_u = N * v;
    timer.stop();
    if (printOutput) {
        std::cout << "Result : " << std::endl;
        printVector(result_col_u);
    }
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Perform norm computations
    // One norm
    timer.start();
    norm = N.norm<algebra::ONE>();
    timer.stop();
    std::cout << "ONE norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Infty norm
    timer.start();
    norm = N.norm<algebra::INFTY>();
    timer.stop();
    std::cout << "INFTY norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Frobenius norm
    timer.start();
    norm = N.norm<algebra::FROBENIUS>();
    timer.stop();
    std::cout << "FROBENIUS norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    std::cout << "=================================================================" << std::endl;

    // Compress the matrix
    N.compress();

    // Perform the multiplication
    std::cout << "Right vector multiplication in compressed state" << std::endl;
    timer.start();
    auto result_col_c = N * v;
    timer.stop();
    if (printOutput) {
        std::cout << "Result : " << std::endl;
        printVector(result_col_c);
    }
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Perform norm computations
    // One norm
    timer.start();
    norm = N.norm<algebra::ONE>();
    timer.stop();
    std::cout << "ONE norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Infty norm
    timer.start();
    norm = N.norm<algebra::INFTY>();
    timer.stop();
    std::cout << "INFTY norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Frobenius norm
    timer.start();
    norm = N.norm<algebra::FROBENIUS>();
    timer.stop();
    std::cout << "FROBENIUS norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    std::cout << "=================================================================" << std::endl;
    
    /********************************************************************************************/

    // Row Major Matrix with complex coefficients
    std::cout << "=================================================================" << std::endl;
    std::cout << "Row Major Matrix with complex coefficients" << std::endl;
    algebra::Matrix<std::complex<double>, algebra::ROW_MAJOR> A("data/test_matrix_complex.mtx");
    if (printOutput) {A.print();}

    // We begin with uncompressed state
    // Perform norm computations
    // One norm
    timer.start();
    norm = A.norm<algebra::ONE>();
    timer.stop();
    std::cout << "ONE norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Infty norm
    timer.start();
    norm = A.norm<algebra::INFTY>();
    timer.stop();
    std::cout << "INFTY norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Frobenius norm
    timer.start();
    norm = A.norm<algebra::FROBENIUS>();
    timer.stop();
    std::cout << "FROBENIUS norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    std::cout << "=================================================================" << std::endl;

    // Compress the matrix
    A.compress();

    // Perform norm computations
    // One norm
    timer.start();
    norm = A.norm<algebra::ONE>();
    timer.stop();
    std::cout << "ONE norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Infty norm
    timer.start();
    norm = A.norm<algebra::INFTY>();
    timer.stop();
    std::cout << "INFTY norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Frobenius norm
    timer.start();
    norm = A.norm<algebra::FROBENIUS>();
    timer.stop();
    std::cout << "FROBENIUS norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    std::cout << "=================================================================" << std::endl;

    /********************************************************************************************/

    // Column Major Matrix with complex coefficients
    std::cout << "=================================================================" << std::endl;
    std::cout << "Column Major Matrix with complex coefficients" << std::endl;
    algebra::Matrix<std::complex<double>, algebra::COL_MAJOR> B("data/test_matrix_complex.mtx");
    if (printOutput) {B.print();} 

    // We begin with uncompressed state
    // Perform norm computations
    // One norm
    timer.start();
    norm = B.norm<algebra::ONE>();
    timer.stop();
    std::cout << "ONE norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Infty norm
    timer.start();
    norm = B.norm<algebra::INFTY>();
    timer.stop();
    std::cout << "INFTY norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Frobenius norm
    timer.start();
    norm = B.norm<algebra::FROBENIUS>();
    timer.stop();
    std::cout << "FROBENIUS norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    std::cout << "=================================================================" << std::endl;

    // Compress the matrix
    B.compress();

    // Perform norm computations
    // One norm
    timer.start();
    norm = B.norm<algebra::ONE>();
    timer.stop();
    std::cout << "ONE norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Infty norm
    timer.start();
    norm = B.norm<algebra::INFTY>();
    timer.stop();
    std::cout << "INFTY norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    // Frobenius norm
    timer.start();
    norm = B.norm<algebra::FROBENIUS>();
    timer.stop();
    std::cout << "FROBENIUS norm : " << norm << std::endl;
    std::cout << "Computation time was " << timer.wallTime() << " microseconds" << std::endl;
    std::cout << "=================================================================" << std::endl;

    return 0;
}