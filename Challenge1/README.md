# 🚀 Challenge 1 – Gradient-Based Minimization of a Multivariate Function
**Authors**: *Marzoli Leo, Murachelli Mattia, Pregnolato Carlo*

---

## 🧠 Project Overview

This project demonstrates the implementation of a solver for gradient-based methods used in the minimization of multivariate functions. The main executable allows the user to choose between several implemented methods and provides the option to use either default solver settings or custom configurations.

---

## 🧮 Available Methods

The following gradient-based minimization methods are implemented:

- **Exponential Decay** [0]
- **Inverse Decay** [1]
- **Armijo** [2]
- **Heavy-Ball** [3]

---

## 🚀 Running the Executable

The executable must be run from the root directory of the project. To launch the program, use the following command:

```bash
./main
```
During execution, the user will be prompted to select the desired minimizer method and choose whether to input custom solver parameters or use the default settings.

## 🛠️ How to Use

A `Makefile` is provided to manage the compilation and cleanup processes efficiently.

### 🧾 Available Make Commands

| Command                | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `make` or `make all`   | Builds the main executable                                                  |
| `make clean`           | Removes all object files and dependencies                                   |
| `make distclean`       | Performs full cleanup, including the removal of compiled executables        |

---

## 📁 Directory Structure

```
.
├── LICENSE
├── Makefile
├── README.md
├── headers
│   ├── Min_arm.hpp
│   ├── Min_exp.hpp
│   ├── Min_hvy.hpp
│   ├── Min_inv.hpp
│   ├── Minimizer.hpp
│   └── vectorExp.hpp
├── main.cpp
└── src
    ├── Min_arm.cpp
    ├── Min_exp.cpp
    ├── Min_hvy.cpp
    ├── Min_inv.cpp
    ├── Minimizer.cpp
    └── vectorExp.cpp
```

---

## 📌 Notes

- All executables are generated in the root project directory.
- The `obj/` folder is automatically created and used to store object (`.o`) and dependency (`.d`) files.