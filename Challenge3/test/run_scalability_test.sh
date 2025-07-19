#!/bin/bash
# SCALABILITY TEST FOR THE LAPLACE SOLVER

# Create a directory to store the results
mkdir -p Results/

# Execution parameters
NUM_TEST=3
REFINEMENT_VECTOR="6 7 8"

# Determine max available physical cores
MAX_CORES=$(($(nproc)/2))

# Test 1 : Sequential
./../main_sequential $NUM_TEST $REFINEMENT_VECTOR

# Test 2 : Parallel 
# Start with 2 cores, increment by 2, up to min(6, MAX_CORES)
for (( CORES=2, COUNT=0; CORES<=6 && COUNT<3 && CORES<=MAX_CORES; CORES+=2, COUNT++ ))
do  
    NUM_THREADS=$((CORES*2))
    OMP_NUM_THREADS=$NUM_THREADS mpirun -np $CORES  ./../main_parallel $NUM_TEST $REFINEMENT_VECTOR
done
