#include "LaplaceSolver.hpp"
#include "Mesh.hpp"

#include "chrono.hpp"
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <omp.h>

using func = LaplaceSolver::func;
using matrix = LaplaceSolver::matrix;

constexpr double pi = 3.14159265358979323846;

int main(int argc, char **argv) {

    // Create parallel environment
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Instantiate OpenMP environment just to retrieve the number of threads
    int num_threads;
    #pragma omp parallel
    {
        #pragma omp master
        {
            num_threads = omp_get_num_threads();
        }
    }

    func f = [](double x, double y) { return 8*pi*pi*std::sin(2*pi*x)*std::sin(2*pi*y); };                  // Forcing term (rhs)
    func g = [](double x, double y) { return 0.0; };                                                        // Dirichlet Boundary Condition
    func uex = [](double x, double y) { return std::sin(2*pi*x)*std::sin(2*pi*y); };                        // Exact Solution

    // Retrieve execution parameters from command line
    int numTests;
    std::vector<std::size_t> refinement_vector;

    numTests = std::stoi(argv[1]);
    if (numTests < 1) {
        std::cout << "Number of tests must be atleast 1" << std::endl;
        return 1;
    }
    refinement_vector.resize(numTests);
    for(int i=0; i<numTests; ++i) {
        refinement_vector[i] = std::stoi(argv[i+2]);
    }
    
    // PARALLEL SOLVER
    // Here rank 0 handles the output, while other processes work together to compute the solution
    if (rank == 0) {
        std::cout << "=====================================================================================" << std::endl;
        std::cout << "======================= Executing on " << size << " processes and " << num_threads << " threads ======================" << std::endl;
        std::cout << "=====================================================================================" << std::endl;
        std::cout << "================================= PARALLEL SOLVER ===================================" << std::endl;
        std::cout << "=====================================================================================" << std::endl;
        }
    // Loop over the refinement vector
    for (std::size_t k=0 ; k < refinement_vector.size(); ++k) {
        // Create the LaplaceSolver
        LaplaceSolver solver(std::pow(2, refinement_vector[k]), f, g, uex);
        // Solve the problem in parallel
        solver("PARALLEL");
        if (rank == 0){
            // Print mesh size
            std::cout << "Mesh size for n = " << std::pow(2, refinement_vector[k]) << ": " << solver.get_mesh_size() << std::endl;
            // Print the execution time
            std::cout << "Execution time for n = " << std::pow(2, refinement_vector[k]) << ": " << solver.get_execution_time() << " seconds" << std::endl;
            // Print the L2 error
            std::cout << "L2 error for n = " << std::pow(2, refinement_vector[k]) << ": " << solver.get_error() << std::endl;
            // Save the solution
            solver.export_to_vtk("PARALLEL");
            std::cout << "-------------------------------------------------------------------------------------" << std::endl;
        }
    }
    if(rank == 0) {
        std::cout << "=====================================================================================" << std::endl;
    }
   
   // Terminate MPI environment
   MPI_Finalize();

   return 0;
}