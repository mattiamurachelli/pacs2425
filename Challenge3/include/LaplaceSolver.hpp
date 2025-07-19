#ifndef LAPLACE_SOLVER_HPP
#define LAPLACE_SOLVER_HPP

#include "Mesh.hpp"

#include <iostream>
#include <functional>
#include <string>

class LaplaceSolver{
    public:
        // Type aliases
        using matrix = std::vector<std::vector<double>>;
        using func = auto(*)(double, double) -> double;
        
        // Constructor
        LaplaceSolver(int n, func f, func g, func uex, std::size_t max_iter = 1e6, double tol = 1e-6);
        // Destructor
        ~LaplaceSolver() = default;
        // Solver
        void operator() (std::string execution_policy);
        // Getters
        double get_error() const;                                                           // Retrieve L2 error
        double get_execution_time() const;                                                  // Retrieve the execution time
        double get_mesh_size() const;                                                       // Retrieve the mesh size
        // Save solution
        void export_to_vtk(std::string) const;                                              // Export the solution to .tvk format (Open in Paraview)

    private:
        // Problem Data
        std::size_t n;                                                                      // Number of rows and columns of the square domain
        Mesh Th;                                                                            // Square mesh (L=1)
        matrix u;                                                                           // Solution
        matrix uold;                                                                        // Previous iteration solution
        func f;                                                                             // Forcing term (rhs)
        func g;                                                                             // Dirichlet Boundary Condition

        // Output evaluation
        func uex;                                                                           // Exact Solution
        double L2error;                                                                     // Error in L2 norm 
        
        // Execution Parameters
        std::size_t max_iter;                                                               // Maximum number of iterations
        double tol;                                                                         // Convergence tolerance
        double execution_time;                                                              // Execution Time

        // Solver functions
        void solve_sequential();                                                            // Sequential Solver
        void solve_parallel();                                                              // Parallel Solver

        // Solver accessory functions
        void apply_BC();                                                                    // Apply Boundary Conditions
        bool check_convergence() const;                                                     // Check Convergence of the method
        bool check_parallel_convergence(std::vector<double>&, std::vector<double>&) const;  // Check Convergence of the parallel method
        void compute_error();                                                               // Compute L2 error
        double compute_local_error(std::vector<double>&, std::vector<double>&) const;       // Compute L2 error of a subregion 
};

#endif // LAPLACE_SOLVER_HPP