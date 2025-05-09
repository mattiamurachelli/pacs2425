# 🚀 Challenge 2 – A sparse Matrix 
**Authors**: *Marzoli Leo, Murachelli Mattia, Pregnolato Carlo*

---

## 🧠 Project Overview

This project showcases the implementation of a templated `Matrix` class in C++, supporting both **row-major** and **column-major** storage orders. The class is generic in terms of the input data type and offers a suite of classical matrix operations, along with:

- A **friend operator** for right vector multiplication
- A templated `norm()` method supporting:
  - **One norm**
  - **Infinity norm**
  - **Frobenius norm**

To evaluate the performance impact of using compressed (`CSC`/`CSR`) vs. uncompressed (`COO`) matrix formats, we use execution timing via the `Chrono.hpp` utility provided by Professor Formaggia.

---

## 🧪 Testing & Demonstration

- **`testing_rowmajor.cpp`** and **`testing_colmajor.cpp`**: Demonstrate basic functionality of the `Matrix` class using both storage formats.
- **`main.cpp`**: Focuses on performance benchmarking of matrix-vector multiplication and norm calculations across different formats and storage orders.

---

## 🧮 Available Matrices
These are the matrix files available in the `data/` folder. Each one is stored in `Matrix Market format (.mtx)`, and can be used as input for testing or benchmarking the implemented matrix operations.
To use any of these matrices, simply pass the filename as a command-line argument when running one of the executables.

- **[lnsp_131](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.html)** : the one assigned in Challenge 2
- **test_matrix_complex** : created by us
- **test_matrix{1/2/3}** : created by us

Furthermore, it is possible to use the executable `MatrixBuilder` (look at next section) to create many more matrices!

---

## 🚀 Running the Executables

All executables must be run from the project root directory. They require command-line arguments to function correctly:

### `main`
Used for performance benchmarking and norm evaluation.

#### Usage :
```bash
./main <matrix_filename> [Y]
```
- **`<matrix_filename>`** : The name of the `.mtx` matrix file in the `data/` folder
- **Y** (optional) : Enables `verbose output`, showing detailed computation steps and intermediate results.
#### Example :
```bash
./main lnsp_131.mtx Y
```

### `test_rowmajor` and `test_colmajor`
Used for functional testing of the Matrix class in row-major and column-major layouts, respectively.
They also offer the possibility of testing matrix-matrix multiplication if two arguments are passed.

#### Usage :
```bash
./test_rowmajor <matrix_filename1> [matrix_filename2]
./test_colmajor <matrix_filename1> [matrix_filename2]
```
- **`<matrix_filename1>`** : The name of the first `.mtx` matrix file in the `data/` folder
- **`<matrix_filename2>`** : The name of the second `.mtx` matrix file in the `data/` folder

#### Example : 
```bash
./test_rowmajor test_matrix.mtx
./test_colmajor test_matrix2.mtx test_matrix3.mtx
```

### `matrixBuilder`
Used for creating random matrices in `.mtx` format.

#### Usage :
```bash
./matrixBuilder <matrix_name> <rows> <cols> <nnz>
```
- **`<matrix_name>`** : The `name of the matrix` that will be stored in MatrixMarket format in the data/ folder
- **`<rows>`** : The `number of rows` of the matrix
- **`<cols>`** : The `number of columns` of the matrix
- **`<nnz>`** : The `number of non-zero-elements` of the matrix

#### Example :
```bash
./matrixBuilder test 100 100 50
```
---

## 🛠️ How to Use

A `Makefile` is provided to manage the compilation and cleanup processes efficiently.

### 🧾 Available Make Commands

| Command                | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `make` or `make all`   | Builds all three executables: `main`, `test_rowmajor`, and `test_colmajor`  |
| `make main_exec`       | Compiles only the `main` executable (performance benchmark)                 |
| `make testing_exec`    | Compiles the two testing executables                                        |
| `make builder_exec`    | Compiles the matrixBuilder executable                                       |
| `make main_clean`      | Cleans object files and dependencies related to the `main` executable       |
| `make testing_clean`   | Cleans object files and dependencies for the testing executables            |
| `make builder_clean`   | Cleans object files and dependencies for the matrixBuilder executable       |
| `make clean`           | Removes all object files and dependencies                                   |
| `make distclean`       | Performs full cleanup, including the removal of compiled executables        |

---

## 📁 Directory Structure

```
.
├── LICENSE
├── Makefile
├── README.md
├── data
│   ├── lnsp_131.mtx
│   ├── test_matrix.mtx
│   ├── test_matrix2.mtx
│   ├── test_matrix3.mtx
│   └── test_matrix_complex.mtx
├── include
│   ├── Matrix.hpp
│   ├── MatrixImplementation.hpp
│   └── chrono.hpp
├── instructions
│   ├── Challenge24-25-2.pdf
│   └── TO_DO.txt
└── src
    ├── main.cpp
    ├── testing_colmajor.cpp
    ├── testing_rowmajor.cpp
    └── RandomMatrixBuilder.cpp
```

---

## 📌 Notes

- All executables are generated in the root project directory.
- The `obj/` folder is automatically created and used to store object (`.o`) and dependency (`.d`) files.
- The `Matrix` class is header-only aside from its main logic and testing files.