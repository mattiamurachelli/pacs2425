// THIS FILE IS USED IN ORDER TO PERFORM TESTS ON THE MATRIX CLASS FOR COL MAJOR MATRIX
#include "Matrix.hpp"
#include "MatrixImplementation.hpp"

int main(int argc, char* argv[]){

    // Check if the user has launched the program correctly
    if (argc == 2){
        std::string filename1 = "data/" + std::string(argv[1]);
        Matrix<double, COL_MAJOR> A(filename1);
        // ATTENTION : Comments on the side are referred to the test with test_matrix.mtx
        // Print the matrix
        A.print();                                                      // Matrix should be printed in uncompressed state

        // Now compress the matrix and print it again
        A.compress();
        A.print();                                                      // Matrix should be printed in compressed state

        // Print the matrix in dense format
        A.printdense();                                                // Matrix should be printed in dense format
        
        // Set a value in the matrix
        A.uncompress();                                                 // Uncompress the matrix
        A.set(1, 3, 6.0);                                               // Adding a value to the matrix

        // Print the matrix again
        A.print();                                                      // Matrix size should also have a new column
        
        // Set another value in the matrix
        A.set(3, 1, 7.0);                                               // Adding a value to the matrix

        // Print the matrix again
        A.print();                                                      // Matrix size should also have a new row

        // Resize the matrix
        A.resize(2, 3);

        // Print the matrix again
        A.print();                                                      // Now the matrix should be resized to 2x3
        
        // Transpose the matrix      
        A.transpose();

        // Print the matrix again
        A.print();                                                      // Now matrix should be transposed

        // Resize, compress and print the matrix again
        A.resize(10,5);
        A.compress();
        A.print();
    }
    else if(argc == 3){
        // Read the two matrices from the files
        std::string filename1 = "data/" + std::string(argv[1]);
        std::string filename2 = "data/" + std::string(argv[2]);
        Matrix<double, COL_MAJOR> A(filename1);
        Matrix<double, COL_MAJOR> B(filename2);

        // Print the matrices
        std::cout << "A: " << std::endl;
        A.print();
        std::cout << "B: " << std::endl;
        B.print();

        // Manipulate the matrices to test the algorithms
        A.resize(5,10);
        B.resize(10,10);

        // Perform matrix multiplication A*B
        Matrix<double, COL_MAJOR> C = A * B;
        std::cout << "C: " << std::endl;
        C.print();

        // Compress the matrices
        A.compress();
        B.compress();

        // Matrix multiplication A*B
        C = A * B;
        C.uncompress();
        C.print();
    }
    else{
        std::cerr << "Usage should be either ./test_colmajor <filename> to test basic class functions or ./test_colmajor <filename1> <filename2> to test matrix product" << std::endl;
        return -1;
    }

    return 0;
}