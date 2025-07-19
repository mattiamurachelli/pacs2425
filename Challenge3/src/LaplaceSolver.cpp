#include "LaplaceSolver.hpp"
#include "chrono.hpp"

#include <math.h>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <omp.h>

// Constructor
LaplaceSolver::LaplaceSolver(int n_, func f_, func g_, func uex_, std::size_t max_iter_, double tol_) :
    n(n_), Th(Mesh(n_)), u(n_, std::vector<double>(n_, 0.0)), uold(n_, std::vector<double>(n_, 0.0)),
    f(f_), g(g_), uex(uex_), L2error(0.0), max_iter(max_iter_), tol(tol_), execution_time(0.0) {}

// Solver
void LaplaceSolver::operator() (std::string execution_policy) {
    // Generate the timer
    Timings::Chrono timer;
    // Apply buondary conditions
    apply_BC();
    // Start the timer 
    timer.start();
    // Solve the problem
    if (execution_policy == "PARALLEL") {
        solve_parallel();
    }
    else if (execution_policy == "SEQUENTIAL") {
        solve_sequential();
    }
    else {
        std::cerr << "Unknown execution policy: " << execution_policy << std::endl;
        return;
    }
    // Stop the timer
    timer.stop();
    // Compute execution time (converted from microsecond to seconds)
    execution_time = 1e-6*timer.wallTime();
    // Compute the error
    compute_error();
}

// Sequential Solver
void LaplaceSolver::solve_sequential() {
    // Iterate until convergence
    for (std::size_t k=0; k<max_iter; ++k) {
        // Update solution
        for (std::size_t i=1; i<n-1; ++i) {
            for (std::size_t j=1 ; j<n-1; ++j) {
                u[i][j] = 0.25 * (uold[i-1][j] + uold[i+1][j] + uold[i][j-1] + uold[i][j+1] + Th.get_h()*Th.get_h()*f(Th(i,j).first, Th(i,j).second));
            }
        }
        /*
        // DEBUG (PRINT ITERATION NUMBER)
        std::cout << "Iteration : " << k << std::endl; 
        // END DEBUG
        */
        // Check convergence
        if (check_convergence()) {
            std::cout << "Converged in " << k << " iterations" << std::endl;
            break;
        }
        // Update old solution for next iteration
        uold = u;
        /*
        // DEBUG (PRINT CURRENT SOLUTION)
        for(std::size_t i = 0; i < n; ++i) {
            for(std::size_t j = 0; j < n; ++j) {
                std::cout << std::setprecision(3) << u[i][j] << "\t";
            }
            std::cout << std::endl;
        }
        // END DEBUG
        */
    }

}

// Parallel Solver
void LaplaceSolver::solve_parallel() {
    // Parallel environment already launched in main.cpp!

    // Define the communicator
    MPI_Comm comm = MPI_COMM_WORLD;

    // Understand our place in the parallel environment
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // We check convergence only every fixed amount of iterations to save some
    // time due to communication overhead
    std::size_t check_conv_iter = 10;

    // Declare on every process the variables that will help us handle
    // communication of the solution between different processes
    std::vector<int> send_idx(size);
    std::vector<int> send_count(size);
    std::vector<int> recv_idx(size);
    std::vector<int> recv_count(size); 

    // Accessory variables for the computation
    // We compute them on every rank since we have less overhead wrt broadcasting
    int remainder = n%size, this_rank_rows, start_idx = 0;
    for(int i=0; i<size; ++i) {
        // Number of rows this process should COMPUTE
        this_rank_rows = (i < remainder) ? (n/size + 1) : (n/size);
        // Sending first row index
        send_idx[i] = start_idx;
        // Number of values to send (ATTENTION! : there are "ghost rows")
        if(i != 0 && i != size-1) {
            send_count[i] = (this_rank_rows + 2)*n;
        } else if(i == 0 || i == size-1) {
            send_count[i] = (this_rank_rows + 1)*n;
        }
        // Receiving first row index
        if(i == 0) {
            recv_idx[i] = 0;
        } else {
            recv_idx[i] = send_idx[i] + n;
        }
        // Number of elements to receive
        recv_count[i] = this_rank_rows*n;
        // Update the starting index for the next process
        if(i == 0) {
            start_idx += n*(this_rank_rows - 1);
        } else {
            start_idx += n*(this_rank_rows);
        }
    }       
    

    // Before scattering our initial solution to all processes, we should flatten
    // it, in order for it to be in contiguous memory slots
    std::vector<double> flattened_global_u;
    if(rank == 0) {
        flattened_global_u.resize(n*n);
        for(std::size_t i = 0; i < n; ++i) {
            for(std::size_t j = 0; j < n; ++j) {
                flattened_global_u[i*n + j] = u[i][j];
            }
        }
    }

    // Declare on each process a vector to store the received data
    std::vector<double> flattened_local_u(send_count[rank]);

    // Declare variables that we will use as support during the computation
    std::size_t num_local_rows = flattened_local_u.size()/n;
    int idx;

    // Now rank 0 can scatter the initial solution to all ranks
    MPI_Scatterv (flattened_global_u.data(), send_count.data(), send_idx.data(), MPI_DOUBLE, flattened_local_u.data(),
        flattened_local_u.size(), MPI_DOUBLE, 0, comm);
    
    // We also declare the old solution for the iterative method updates
    std::vector<double> flattened_local_uold = flattened_local_u;

    // And also the MPI_Requests for thr ghost rows exchange with non-blocking communications
    MPI_Request send_requests[2], recv_requests[2];
    
    // We also need to remember that, even though every process has access to the full mesh, it has to work on
    // its assigned subregion. We therefore compute the row offset of each process wrt the first mesh row
    int mesh_row_offset = send_idx[rank]/n;

    // Now we should start computing the solution with the iterative method
    for(std::size_t it=0; it < max_iter; ++it) {
        
        // If not converged we need to exchange values for ghost rows before proceding
        // with the next iteration. We use non-blocking communications to be able
        // to overlap the communication with the computation of the interior values
        // that do not require to wait for the ghost rows exchange

        // All ranks except the first one pass their first non-ghost row as the last (ghost)
        // row of the previous rank, while all ranks except the last one pass their last
        // non-ghost row as the first (ghost) row of the next rank
        /*
        if(rank > 0) { 
            //MPI_Sendrecv (&flattened_local_uold[n], n, MPI_DOUBLE, rank-1, 0, flattened_local_uold.data(),
            //                n, MPI_DOUBLE, rank-1, 1, comm, MPI_STATUS_IGNORE);
            MPI_Isend (&flattened_local_uold[n], n, MPI_DOUBLE, rank-1, 0, comm, &send_requests[0]);
            MPI_Irecv (&flattened_local_uold[(num_local_rows-1)*n], n, MPI_DOUBLE, rank+1, 1, comm, &recv_requests[0]);

        }
        if(rank < size - 1){
            //MPI_Sendrecv (&flattened_local_uold[(num_local_rows-2)*n], n, MPI_DOUBLE, rank+1, 1, &flattened_local_uold[(num_local_rows-1)*n],
            //                n, MPI_DOUBLE, rank+1, 0, comm, MPI_STATUS_IGNORE);
            MPI_Isend (&flattened_local_uold[(num_local_rows-2)*n], n, MPI_DOUBLE, rank-1, 1, comm, &send_requests[1]);
            MPI_Irecv (flattened_local_uold.data(), n, MPI_DOUBLE, rank-1, 0, comm, &recv_requests[1]);
        }
        */
       if(rank > 0) { 
            MPI_Isend(&flattened_local_uold[n], n, MPI_DOUBLE, rank-1, 0, comm, &send_requests[0]);
            MPI_Irecv(flattened_local_uold.data(), n, MPI_DOUBLE, rank-1, 1, comm, &recv_requests[0]);
        }
        if(rank < size - 1){
            MPI_Isend(&flattened_local_uold[(num_local_rows-2)*n], n, MPI_DOUBLE, rank+1, 1, comm, &send_requests[1]);
            MPI_Irecv(&flattened_local_uold[(num_local_rows-1)*n], n, MPI_DOUBLE, rank+1, 0, comm, &recv_requests[1]);
        }

        // Compute the values for the local solution using OpenMP multithreading only in the interior of the subregion
        // while we wait for the ghost rows exchange to be completed
        #pragma omp parallel for collapse(2) private(idx)
        for(std::size_t i = 2; i < num_local_rows - 2; ++i) { 
            for(std::size_t j = 1 ; j < n - 1; ++j) {
                idx = i*n + j;
                flattened_local_u[idx] = 0.25*(flattened_local_uold[idx-n] + flattened_local_uold[idx+n]
                                       + flattened_local_uold[idx-1] + flattened_local_uold[idx+1]
                                       + Th.get_h()*Th.get_h()*f(Th(i+mesh_row_offset,j).first, Th(i+mesh_row_offset,j).second));
            }
        }

        // Now we need to make sure that the ghost rows exchange has been properly executed
        if(rank > 0){
            MPI_Wait(&send_requests[0], MPI_STATUS_IGNORE);
            MPI_Wait(&recv_requests[0], MPI_STATUS_IGNORE);
        }
        if(rank < size - 1){
            MPI_Wait(&send_requests[1], MPI_STATUS_IGNORE);
            MPI_Wait(&recv_requests[1], MPI_STATUS_IGNORE);
        }

        // Now we can also compute the boundary rows
        #pragma omp parallel for private(idx)
            for(std::size_t j = 1; j < n - 1; ++j){
                std::size_t i = 1;
                idx = i*n + j;
                flattened_local_u[idx] = 0.25*(flattened_local_uold[idx-n] + flattened_local_uold[idx+n]
                                       + flattened_local_uold[idx-1] + flattened_local_uold[idx+1]
                                       + Th.get_h()*Th.get_h()*f(Th(i+mesh_row_offset,j).first, Th(i+mesh_row_offset,j).second));
            }

        #pragma omp parallel for private(idx)
            for(std::size_t j = 1; j < n - 1; ++j){
                std::size_t i = num_local_rows-2;
                idx = i*n + j;
                flattened_local_u[idx] = 0.25*(flattened_local_uold[idx-n] + flattened_local_uold[idx+n]
                                       + flattened_local_uold[idx-1] + flattened_local_uold[idx+1]
                                       + Th.get_h()*Th.get_h()*f(Th(i+mesh_row_offset,j).first, Th(i+mesh_row_offset,j).second));
            }

        // Perform a local convergence test to see if we found the solution already
        if(it % check_conv_iter == 0 && check_parallel_convergence(flattened_local_u, flattened_local_uold)) {
            if(rank == 0) {
                std::cout << "Converged in " << it << " iterations" << std::endl;
            }
            break;
        }
        /*
        // DEBUG (PRINT ITERATION NUMBER)
        else {
            if(rank == 0){
                std::cout << "Iteration : " << it << std::endl;
            }
        }
        //END DEBUG
        */
       /*
        // DEBUG (PRINT CURRENT SOLUTION FROM DIFFERENT RANKS)
        for(int i = 0; i < size; ++i) {
            MPI_Barrier(comm);
        
            if(rank == i) {
                std::cout << "Rank " << rank << ":\n";
                int count = 0;
                for(std::size_t k = 0; k < flattened_local_u.size(); ++k) {
                    std::cout << std::setprecision(3) <<flattened_local_u[k] << "\t";
                    count++;
                    if(count % n == 0) std::cout << "\n";
                }
                std::cout << std::endl;
            }
        
            MPI_Barrier(comm); 
        }
        // END DEBUG
        */
        /*
        // DEBUG (PRINT CURRENT SOLUTION)
        // Assemble the final solution of the problem
        MPI_Gatherv ((rank == 0) ? flattened_local_u.data() : &flattened_local_u[n], recv_count[rank], MPI_DOUBLE, flattened_global_u.data(),
                    recv_count.data(), recv_idx.data(), MPI_DOUBLE, 0, comm);

        // De-flatten it on rank 0
        if(rank == 0) {
            for(std::size_t i = 0; i < n; ++i) {
                for(std::size_t j = 0; j < n; ++j) {
                    u[i][j] = flattened_global_u[i*n + j];
                }
            }
            
            std::cout << "Current global solution : " << std::endl;

            for(std::size_t i = 0; i < n; ++i) {
                for(std::size_t j = 0; j < n; ++j) {
                    std::cout << std::setprecision(3) << u[i][j] << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        // END DEBUG
        */
        // Update old solution for next iteration
        std::swap(flattened_local_uold, flattened_local_u);

    }
    
    // Assemble the final solution of the problem
    MPI_Gatherv ((rank == 0) ? flattened_local_u.data() : &flattened_local_u[n], recv_count[rank], MPI_DOUBLE, flattened_global_u.data(),
                recv_count.data(), recv_idx.data(), MPI_DOUBLE, 0, comm);

    // De-flatten it on rank 0
    if(rank == 0) {
        for(std::size_t i = 0; i < n; ++i) {
            for(std::size_t j = 0; j < n; ++j) {
                u[i][j] = flattened_global_u[i*n + j];
            }
        }
        /* 
        // DEBUG (PRINT FINAL SOLUTION)
        std::cout << "Final solution : " << std::endl;
        for(std::size_t i = 0; i < n; ++i) {
            for(std::size_t j = 0; j < n; ++j) {
                std::cout << std::setprecision(3) << u[i][j] << "\t";
            }
            std::cout << std::endl;
        }
        //END DEBUG
        */
    }
}

// Apply Boundary Conditions
void LaplaceSolver::apply_BC() {
    // Top side
    for (std::size_t j=0; j<n; ++j) {
        u[0][j] = g(Th(0,j).first, Th(0,j).second);
    }
    // Left side (excluding corners)
    for (std::size_t i=1; i<n-1; ++i) {
        u[i][0] = g(Th(i,0).first, Th(i,0).second);
    }
    // Bottom side
    for (std::size_t j=0; j<n; ++j) {
        u[n-1][j] = g(Th(n-1,j).first, Th(n-1,j).second);
    }
    // Right side (excluding corners)
    for (std::size_t i=1; i<n-1; ++i) {
        u[i][n-1] = g(Th(i,n-1).first, Th(i,n-1).second);
    }
    // Apply the conditions to uold too
    uold = u;
}

// Check Convergence of the method
bool LaplaceSolver::check_convergence() const {
    // Compute error
    double e = 0.0;
    for (std::size_t i=0; i<n; ++i){
        for (std::size_t j=1; j<n-1; ++j){
            e += (u[i][j] - uold[i][j])*(u[i][j] - uold[i][j]);
        }
    }
    e = std::sqrt(Th.get_h() * e);
    /*
    // DEBUG (PRINT ERROR)
    std::cout << "Error : " << e << std::endl;
    // END DEBUG
    */
    // Check convergence wrt tolerance
    return (e < tol) ? true : false;
}

// Check convergence of the parallel method
bool LaplaceSolver::check_parallel_convergence(std::vector<double> &u, std::vector<double> &uold) const {
    // Compute local error
    double local_error = compute_local_error(u, uold);

    // Retrieve global error
    double global_error;
    MPI_Allreduce (&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    global_error = std::sqrt(Th.get_h()*global_error);
    /*
    // DEBUG (PRINT GLOBAL ERROR)
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0){
        std::cout << "Global error =  " << global_error << std::endl;
    }
    //END DEBUG
    */
    // Check wrt tolerance value
    return (global_error < tol) ? true : false;
}

// Compute error of a subregion for the parallel method
double LaplaceSolver::compute_local_error(std::vector<double> &u, std::vector<double> &uold) const {
    // Compute error
    double e = 0.0;
    std::size_t num_local_rows = u.size()/n;
    int idx;
    // We exclude first and last row since for every process they are 
    // either boundary rows or ghost rows
    for(std::size_t i = 1; i < num_local_rows - 1; ++i) {
        for(std::size_t j = 1; j < n-1 ; ++j) {
            idx = i*n + j;
            e += (u[idx] - uold[idx])*(u[idx] - uold[idx]);
        }
    }
     return e;
}

// Compute L2 error
void LaplaceSolver::compute_error() {
    for(std::size_t i=0; i<n; ++i) {
        for(std::size_t j=0; j<n; j++) {
            L2error += (u[i][j] - uex(Th(i,j).first, Th(i,j).second))*(u[i][j] - uex(Th(i,j).first, Th(i,j).second));
        }
    }
    L2error = std::sqrt(Th.get_h() * L2error);
}

// Function to export the solution to .vtk format
void LaplaceSolver::export_to_vtk(std::string execution_policy) const {
    // Open the file
    std::string filename;
    if(execution_policy == "SEQUENTIAL") {
        filename = "Results/solution_" + std::to_string(n) + "_sequential.vtk";
    } else if(execution_policy == "PARALLEL"){
        filename = "Results/solution_" + std::to_string(n) + "_parallel.vtk";
    } else {
        std::cout << "Unknown execution policy, solution wasn't generated properly!" << std::endl;
        return;
    }
    std::ofstream file(filename);
    // Check if the file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing" << std::endl;
        return;
    }
    // Write the VTK header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "Laplace Solver Solution" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET STRUCTURED_POINTS" << std::endl;
    file << "DIMENSIONS " << n << " " << n << " 1" << std::endl;
    file << "ORIGIN 0 0 0" << std::endl;
    file << "SPACING " << Th.get_h() << " " << Th.get_h() << " 1" << std::endl;
    file << "POINT_DATA " << n*n << std::endl;
    file << "SCALARS solution double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    // Write the solution data (in a row-major order)
    for (std::size_t i=0; i<n; ++i) {
        for (std::size_t j=0; j<n; ++j) {
            file << u[i][j] << std::endl;
        }
    }
    // Close the file
    file.close();
    // Notify the user
    std::cout << "Solution exported to " << filename << std::endl;
}

// L2 Error getter
double LaplaceSolver::get_error() const {
    return L2error;
}

// Execution time getter
double LaplaceSolver::get_execution_time() const {
    return execution_time;
}

// Mesh size getter
double LaplaceSolver::get_mesh_size() const {
    return Th.get_h();
}