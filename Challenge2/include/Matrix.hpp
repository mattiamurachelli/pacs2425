#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <vector>
#include <string>

namespace algebra{

    constexpr size_t MAX_SIZE_PRINTABLE = 50;

    enum StorageOrder
    {
        ROW_MAJOR,
        COL_MAJOR
    };

    enum NormType
    {
        ONE,
        INFTY,
        FROBENIUS
    };

    template<typename T>
    struct CompressedMatrix{
        std::vector<std::size_t> row_index;
        std::vector<std::size_t> col_index;
        std::vector<T> values;
    };

    template<StorageOrder Order>
    struct MatrixOrdering{
        bool operator() (const std::pair<std::size_t, std::size_t> &p1, const std::pair<std::size_t, std::size_t> &p2) const {
            if constexpr (Order == ROW_MAJOR) {
                return p1 < p2; // use std::less
            } else {
                return (p1.second < p2.second) || (p1.second == p2.second && p1.first < p2.first);
            }
        }
    };

    // Forward declarations
    template<typename T, StorageOrder Order>
    class Matrix;

    template<typename T, StorageOrder Order>
    std::vector<T> operator*(const Matrix<T, Order>&, const std::vector<T>&);

    template<typename T, StorageOrder Order>
    Matrix<T,Order> operator*(const Matrix<T, Order>&, const Matrix<T, Order>&);

    // Class definition
    template<typename T, StorageOrder Order>
    class Matrix{
        public :
            // Friend operators
            friend std::vector<T> operator*<>(const Matrix<T, Order>&, const std::vector<T>&);
            friend Matrix<T,Order> operator*<>(const Matrix<T, Order>&, const Matrix<T, Order>&);

            // Type aliases
            using CompressedMatrix = algebra::CompressedMatrix<T>;
            using UncompressedMatrix = std::map<std::pair<std::size_t, std::size_t>, T, MatrixOrdering<Order>>;

            // Constructor and Destructor
            Matrix(std::string);
            Matrix() = default;
            ~ Matrix() = default;

            // Functions
            void set(std::size_t, std::size_t, T);
            void setRows(std::size_t);
            void setCols(std::size_t);
            void setNnz(std::size_t);
            std::size_t getRows() const;
            std::size_t getCols() const;
            std::size_t getNnz() const;
            T operator()(std::size_t, std::size_t) const;
            T& operator()(std::size_t, std::size_t);
            void resize(std::size_t, std::size_t);
            void transpose();
            bool BuyFromMarket(std::string);
            bool SellToMarket(std::string);
            bool is_compressed() const;
            void compress();
            void uncompress();
            void print() const;
            void printdense() const;

            template<NormType N>
            double norm() const;

            private:
            // Data
            UncompressedMatrix uData;
            CompressedMatrix cData;
            std::size_t rows = 0;
            std::size_t cols = 0;
            std::size_t nnz = 0;
            bool compressed = 0;
            
    };

    template<typename T, StorageOrder Order>
    std::vector<T> operator*(const Matrix<T, Order>&, const std::vector<T>&);

    template<typename T, StorageOrder Order>
    Matrix<T,Order> operator*(const Matrix<T, Order>&, const Matrix<T, Order>&);
}

#endif // MATRIX_HPP