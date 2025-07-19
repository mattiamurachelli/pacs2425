#!/bin/bash
# SCALABILITY TEST FOR THE LAPLACE SOLVER

# Create a directory to store the results
mkdir -p Results/

# Execution parameters
NUM_TEST=3
REFINEMENT_VECTOR="7 8 9"

# Determine max available physical cores
MAX_CORES=$(($(nproc)/2))

# Test 1 : Sequential
./../main_sequential $NUM_TEST $REFINEMENT_VECTOR

# Test 2 : Parallel 
# Use only two cores and the maximum number of threads allowed
OMP_NUM_THREADS=$(nproc) mpirun -np 2 ./../main_parallel $NUM_TEST $REFINEMENT_VECTOR