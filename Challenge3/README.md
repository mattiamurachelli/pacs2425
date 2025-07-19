# ⚙️ Challenge 3 – A matrix-free parallel solver for the Laplace equation
**Authors**: *Marzoli Leo, Murachelli Mattia, Pregnolato Carlo*

---

## 🧠 Project Overview

In this challenge, we implemented a solver for the Laplace equation with Dirichlet boundary conditions. The solver supports both:

- Sequential execution
- Parallel execution (via `MPI` and `OpenMP`)

The program solves the equation over increasingly refined structured meshes, showcasing scalability with a variable number of processors.
This project allows for performance benchmarking and direct visualization of the computed solutions. A dedicated bash script automates the scalability tests.

---

## 🚀 Running the Solver

### 🔨 Compilation
From the project root directory, compile everything using the provided `Makefile`:
```bash
make
```
This will generate both sequential and parallel executables.

### 📁 Test & Scalability Script
Move into the `test/` directory and run the test script:
```bash
cd test/
./run_scalability_test.sh
```
This script:
- Solves the Laplace equation sequentially
- Then solves it in parallel, scaling the number of processors
- Uses several mesh refinements
- Stores output solutions in the Results/ folder

🛠️ Note: If the script is not executable, make it so with:
```bash
chmod +x run_scalability_test.sh
```
You are also free to modify the script in order to change the `number of tests` to perform and the `mesh refinement` values.

---

## 📂 Directory Structure
```bash
.
├── Challenge24-25-3.pdf
├── LICENSE
├── Makefile
├── README.md
├── hw.info
├── include
│   ├── LaplaceSolver.hpp
│   ├── Mesh.hpp
│   └── chrono.hpp
├── src
│   ├── LaplaceSolver.cpp
│   ├── Mesh.cpp
│   ├── main_parallel.cpp
│   └── main_sequential.cpp
└── test
    ├── run_scalability_test.sh
    └── run_scalability_test_optimal.sh
```

---

## ⚙️ Makefile Targets

| Command                | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `make` or `make all`   | Builds both the sequential and parallel solvers                             |
| `make clean`           | Cleans all object and intermediate files                                    |
| `make distclean`       | Performs full cleanup, including executables                                |

---

## 📈 Output & Visualization

- Output files are saved in `test/Results/`
- They are `compatible with` standard VTK viewers such as `ParaView`
- Each output corresponds to a specific `refinement level` and `execution policy`

--- 

## 📝 Notes

- All executables are placed in the root project directory
- The `run_scalability_test.sh` script is the easiest way to reproduce and test the scalability
- Ensure that MPI is installed to run the parallel executable (mpirun is required)

## 🧪 Testing Phase

For this final challenge, we decided to create a dedicated testing section, as the results we obtained were, in a sense, `suboptimal` compared to our expectations. Below, we provide a brief analysis of our findings.

Regarding the use of `OpenMP`, the parallelization of the computation using the `parallel for collapse` pragma led to a decent speed-up.

On the other hand, the use of `MPI` resulted in a speed-up only when the `number of processes` was equal to `2`. We suspect this is due to the overhead introduced by ghost-row communication between processes, which becomes impractical when the number of processes is too high.

That said, the `run_scalability_test.sh` script may not be ideal for showcasing the capabilities of our code, as it tends to produce poor performance results. To address this, we created an alternative script, `run_scalability_test_optimal.sh`, which uses only `two processes` and the `maximum number of threads` to demonstrate our program's performance in its most optimized environment.

The results below were obtained using this optimized version on our hardware setup, as described in the `hw.info` file:
```bash
=====================================================================================
================================ SEQUENTIAL SOLVER ==================================
=====================================================================================
Converged in 7220 iterations
Mesh size for n = 128: 0.00787402
Execution time for n = 128: 4.42541 seconds
L2 error for n = 128: 0.000334006
Solution exported to Results/solution_128_sequential.vtk
-------------------------------------------------------------------------------------
Converged in 25670 iterations
Mesh size for n = 256: 0.00392157
Execution time for n = 256: 80.9656 seconds
L2 error for n = 256: 0.00288879
Solution exported to Results/solution_256_sequential.vtk
-------------------------------------------------------------------------------------
Converged in 89296 iterations
Mesh size for n = 512: 0.00195695
Execution time for n = 512: 1219.48 seconds
L2 error for n = 512: 0.0130852
Solution exported to Results/solution_512_sequential.vtk
-------------------------------------------------------------------------------------
=====================================================================================
=====================================================================================
======================= Executing on 2 processes and 12 threads ======================
=====================================================================================
================================= PARALLEL SOLVER ===================================
=====================================================================================
Converged in 7220 iterations
Mesh size for n = 128: 0.00787402
Execution time for n = 128: 4.19873 seconds
L2 error for n = 128: 0.000334006
Solution exported to Results/solution_128_parallel.vtk
-------------------------------------------------------------------------------------
Converged in 25670 iterations
Mesh size for n = 256: 0.00392157
Execution time for n = 256: 31.2298 seconds
L2 error for n = 256: 0.00288879
Solution exported to Results/solution_256_parallel.vtk
-------------------------------------------------------------------------------------
Converged in 89300 iterations
Mesh size for n = 512: 0.00195695
Execution time for n = 512: 385.588 seconds
L2 error for n = 512: 0.0130812
Solution exported to Results/solution_512_parallel.vtk
-------------------------------------------------------------------------------------
=====================================================================================
```

---