#include "Matrix.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <cmath>
#include <complex>
#include <type_traits>
#include <algorithm>

using namespace algebra;

// CONSTRUCTOR
template <typename T, StorageOrder Order>
Matrix<T, Order>::Matrix(std::string filename){
    // If MatrixMarket file is not found, use default constructor
    if(!BuyFromMarket(filename)){
        Matrix();
    }
}

// METHOD USED TO SET VALUES
template <typename T, StorageOrder Order>
void Matrix<T, Order>::set(std::size_t row_index, std::size_t col_index, T value){
    // Matrix must be in uncompressed state to set a value
    bool wasCompressed = false;
    if(compressed){
        this->uncompress();
        wasCompressed = true;
    }
    // we check one dimension at a time
    if(row_index >= rows){
        this->resize(row_index+1, cols);
    }
    if(col_index >= cols){
        this->resize(rows, col_index+1);
    }
    // Now we can set the value (if it is already in the map it will be overwritten)
    // we increase nnz only if the value was not already there
    if (uData.find(std::make_pair(row_index, col_index)) == uData.end()) {nnz++;}
    uData[std::make_pair(row_index, col_index)] = value;
    // Restore the original state of the matrix if it's necessary
    if (wasCompressed) {
        this->compress();
    }
}

// METHOD TO SET NUMBER OF ROWS
template <typename T, StorageOrder Order>
void Matrix<T, Order>::setRows(std::size_t num_rows) {
    rows = num_rows;
}

// METHOD TO SET NUMBER OF COLUMNS
template <typename T, StorageOrder Order>
void Matrix<T, Order>::setCols(std::size_t num_cols) {
    cols = num_cols;
}

// METHOD TO SET NUMBER OF NON ZERO ELEMENTS
template <typename T, StorageOrder Order>
void Matrix<T, Order>::setNnz(std::size_t nnz_) {
    nnz = nnz_;
}

// METHOD USED TO RETRIEVE NUMBER OF ROWS
template <typename T, StorageOrder Order>
std::size_t Matrix<T, Order>::getRows() const{
    return rows;
}

// METHOD USED TO RETRIEVE NUMBER OF COLUMNS
template <typename T, StorageOrder Order>
std::size_t Matrix<T, Order>::getCols() const{
    return cols;
}

// METHOD USED TO RETRIEVE NUMBER OF NON ZERO ELEMENTS
template <typename T, StorageOrder Order>
std::size_t Matrix<T, Order>::getNnz() const{
    return nnz;
}

// METHOD USED TO ACCESS VALUES FROM A CONST MATRIX
template <typename T, StorageOrder Order>
T Matrix<T, Order>::operator()(std::size_t row_index, std::size_t col_index) const{
    // First we check if the indices are out of bounds, then if they are not we look for
    // and return the value 
    if( row_index >= rows || col_index >= cols){
        std::cerr << "Error: Index out of bounds" << std::endl;
        return std::numeric_limits<T>::quiet_NaN();
    }

    if(compressed){
        if(Order == ROW_MAJOR){
            std::size_t start = cData.row_index[row_index];
            std::size_t end   = cData.row_index[row_index+1];
            for(std::size_t i=start; i<end; ++i){
                if(cData.col_index[i] == col_index){
                    return cData.values[i];
                }
            }
            return 0;
        }
        if(Order == COL_MAJOR){
            std::size_t start = cData.col_index[col_index];
            std::size_t end   = cData.col_index[col_index+1];
            for(std::size_t i=start; i<end; ++i){
                if(cData.row_index[i] == row_index){
                    return cData.values[i];
                }
            }
            return 0;
        }
    }
    // The StorageOrder is irrelevant in this case
    else{
        auto it = uData.find(std::make_pair(row_index, col_index));
        // If the value is found, return it, otherwise return 0
        if(it != uData.end()){
            return it->second;
        } else {
            return 0;
        }
    }
}

// METHOD TO RETRIEVE VALUES FROM A MATRIX
// THIS CAN BE MADE ONLY WHEN THE FUNCTION IS UNCOMPRESSED AND IT DOES NOT LET THE
// USER ACCESS VALUES THAT HAVEN'T BEEN SET YET
template <typename T, StorageOrder Order>
T& Matrix<T, Order>::operator()(std::size_t row_index, std::size_t col_index){
    // First we check if the indices are out of bounds, then if they are not we look for
    // and return the value 
    if( row_index >= rows || col_index >= cols){
        std::cerr << "Error: Index out of bounds" << std::endl;
        return std::numeric_limits<T>::quiet_NaN();
    }

    if(compressed){
        if(Order == ROW_MAJOR){
            std::size_t start = cData.row_index[row_index];
            std::size_t end   = cData.row_index[row_index+1];
            for(std::size_t i=start; i<end; ++i){
                if(cData.col_index[i] == col_index){
                    return cData.values[i];
                }
            }
            return 0;
        }
        if(Order == COL_MAJOR){
            std::size_t start = cData.col_index[col_index];
            std::size_t end   = cData.col_index[col_index+1];
            for(std::size_t i=start; i<end; ++i){
                if(cData.row_index[i] == row_index){
                    return cData.values[i];
                }
            }
            return 0;
        }
    }
    // The StorageOrder is irrelevant in this case
    else{
        auto it = uData.find(std::make_pair(row_index, col_index));
        // If the value is found, return it, otherwise return 0
        if(it != uData.end()){
            return it->second;
        } else {
            return 0;
        }
    }
}

// METHOD USED TO RESIZE THE MATRIX
template <typename T, StorageOrder Order>
void Matrix<T, Order>::resize(std::size_t num_rows, std::size_t num_cols){
    // ATTENTION : num_rows and num_cols are the new dimensions of the matrix,
    // they are not intended to be used to access the matrix values (as they are one past the last valid index)

    // Matrix must be in uncompressed state to resize
    bool wasCompressed = false;
    if(compressed){
        this->uncompress();
        wasCompressed = true;
    }
    // If the matrix enlarges there is no problem, we just need to change sizes
    // Othwerwise we need to get rid of the values that are out of bounds
    if(num_rows > rows && num_cols > cols){
        rows = num_rows;
        cols = num_cols;
    } else {
        auto it = uData.cbegin();
        while(it != uData.cend()){
            if(it->first.first >= num_rows || it->first.second >= num_cols){
                it = uData.erase(it);
                nnz--;
            }
            else{
                ++it;
            }
        }
        // Now we can resize the matrix
        rows = num_rows;
        cols = num_cols;
    }
    // Restore the original state of the matrix if it's necessary
    if (wasCompressed) {
        this->compress();
    }
}

// FUNCTION USED TO TRANSPOSE THE MATRIX
template <typename T, StorageOrder Order>
void Matrix<T, Order>::transpose(){
    // The StorageOrder is irrelevant in this case, we simply need to swap the indices
    // As usual, to treat matrix values, we uncompress the matrix first
    bool wasCompressed = false;
    if(compressed){
        this->uncompress();
        wasCompressed = true;
    }
    // Create new matrix
    UncompressedMatrix transposedMatrix;
    // Iterate over the data and insert values into transposed Matrix with swapped indices
    for(auto& [key, value] : uData){
        transposedMatrix[std::make_pair(key.second, key.first)] = value;
    }
    // transpose the matrix
    uData = transposedMatrix;
    std::swap(rows, cols);
    // Restore the original state of the matrix if it's necessary
    if (wasCompressed) {
        this->compress();
    }
}

// METHOD USED TO IMPORT MATRIX VALUES FROM FILE
// IT BASICALLY IS THE CONSTRUCTOR OF THE MATRIX
template <typename T, StorageOrder Order>
bool Matrix<T, Order>::BuyFromMarket(std::string filename){

    std::cout << "Buying the matrix from MatrixMarket..." << std::endl;
    std::ifstream file(filename);
    
    // Check if the file was opened correctly
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        // Return false to indicate failure
        return false;
    }

    std::string line;
    std::size_t nrows, ncols, nonZeros;
    
    // Skip first rows with comments
    while (getline(file, line)) {
        if (line[0] != '%') break;
    }
    
    // Read matrix dimensions and number of non-zero-values from the first uncommented line
    std::istringstream reader(line);
    reader >> nrows >> ncols >> nonZeros;

    // Initialize matrix dimensions and state
    rows = nrows;
    cols = ncols;
    nnz = nonZeros;
    compressed = false;

    // Read the matrix values and fill uData based on the StorageOrder after making sure that the map is empty
    uData.clear();

    std::size_t currRow, currCol;
    T currVal;
    double real, imag;

    if constexpr (std::is_same_v<T, std::complex<double>>) {
        // import the uncompressed data
        while(getline(file,line)){
            std::istringstream reader(line);
            reader >> currRow >> currCol >> real >> imag;
            uData.insert(std::make_pair(std::make_pair(currRow - 1, currCol - 1), std::complex<double>(real, imag)));
        }
    } else {
        // import the uncompressed data
        while(getline(file,line)){
            std::istringstream reader(line);
            reader >> currRow >> currCol >> currVal;
            uData.insert(std::make_pair(std::make_pair(currRow - 1, currCol - 1), currVal));
        }
    }
    
    // Close the file
    file.close();

    std::cout << "Matrix bought from MatrixMarket!" << std::endl;
    
    // Return true to indicate success
    return true;
}

// METHOD USED TO EXPORT A MATRIX TO FILE IN .MTX FORMAT
template <typename T, StorageOrder Order>
bool Matrix<T, Order>::SellToMarket(std::string filename) {

    // Create output file stream
    std::ofstream output("data/" + filename + ".mtx");
    if (!output) {
        std::cerr << "Error: Cannot open file 'data/" << filename << "'\n";
        return false;
    }

    // Matrix will be exported in compressed format
    bool wasCompressed = false;
    if (compressed){
        this->uncompress();
        wasCompressed = true; 
    }

    // Fill the file
    // Header
    output << "%%MatrixMarket matrix coordinate real general\n";
    // Rows Columns NonZeroElements
    output << rows << " " << cols << " " << nnz << "\n";
    // Values
    for (auto & [key, value] : uData) {
        output << key.first+1 << " " << key.second+1 << " " << value << "\n"; 
    }

    // Restore original state if matrix was compressed
    if (wasCompressed) {
        this->compress();
    }

    // Close the file
    output.close();

    std::cout << "Matrix sold to MatrixMarket!" << std::endl;

    // Return true to indicate success
    return true;
}

// METHOD USED TO CHECK IF THE MATRIX IS COMPRESSED
template <typename T, StorageOrder Order>
bool Matrix<T, Order>::is_compressed() const{
    return compressed;
}

// METHOD USED TO COMPRESS THE MATRIX FROM UNCOMPRESSED STATE
template <typename T, StorageOrder Order>
void Matrix<T, Order>::compress(){
    if(is_compressed()){
        std::cout << "Matrix is already compressed." << std::endl;
        return;
    }
    compressed = true;
    if constexpr(Order == ROW_MAJOR){
        cData.row_index.clear();
        cData.col_index.clear();
        cData.values.clear();

        cData.row_index.reserve(rows+1);
        cData.col_index.reserve(uData.size());
        cData.values.reserve(uData.size());

        std::size_t currentRow = 0;
        std::size_t counter = 0;
        cData.row_index.push_back(currentRow);
        for(auto a=uData.cbegin(); a!=uData.cend(); ++a){
            cData.values.push_back(a->second);
            cData.col_index.push_back(a->first.second);
            if(a->first.first == currentRow){
                ++counter;
            }
            else{
                cData.row_index.push_back(counter);
                for(std::size_t i=currentRow+1; i<a->first.first; ++i){
                    cData.row_index.push_back(counter);
                }
                currentRow = a->first.first;
                ++counter;
            }
        }
        while(currentRow < rows){
            cData.row_index.push_back(counter);
            ++currentRow;
        }
    }
    if constexpr(Order == COL_MAJOR){
        cData.row_index.clear();
        cData.col_index.clear();
        cData.values.clear();

        cData.row_index.reserve(uData.size());
        cData.col_index.reserve(cols+1);
        cData.values.reserve(uData.size());

        std::size_t currentCol = 0;
        std::size_t counter = 0;
        cData.col_index.push_back(currentCol);
        for(auto a=uData.cbegin(); a!=uData.cend(); ++a){
            cData.values.push_back(a->second);
            cData.row_index.push_back(a->first.first);
            if(a->first.second == currentCol){
                ++counter;
            }
            else{
                cData.col_index.push_back(counter);
                for(std::size_t i=currentCol+1; i<a->first.second; ++i){
                    cData.col_index.push_back(counter);
                }
                currentCol = a->first.second;
                ++counter;
            }
        }
        while(currentCol < cols){
            cData.col_index.push_back(counter);
            ++currentCol;
        }
    }
    uData.clear();
}

// METHOD USED TO UNCOMPRESS THE MATRIX FROM COMPRESSED STATE
template<typename T, StorageOrder Order>
void Matrix<T, Order>::uncompress() {
    if (!compressed) {
        return; // Already uncompressed
    }

    uData.clear(); // Clear the uncompressed data

    if constexpr(Order == ROW_MAJOR) {
        // Iterate over the rows
        for (std::size_t row = 0; row < cData.row_index.size() - 1; ++row) {
            std::size_t start = cData.row_index[row];
            std::size_t end = cData.row_index[row + 1];

            // Iterate over the non-zero elements in the current row
            for (std::size_t idx = start; idx < end; ++idx) {
                std::size_t col = cData.col_index[idx]; // Column index
                T value = cData.values[idx];           // Non-zero value

                // Insert into COOmap (uData)
                uData[{row, col}] = value;
            }
        }
    } else {
        // For column-major order, iterate over the columns
        for (std::size_t col = 0; col < cData.col_index.size() - 1; ++col) {
            std::size_t start = cData.col_index[col];
            std::size_t end = cData.col_index[col + 1];

            // Iterate over the non-zero elements in the current column
            for (std::size_t idx = start; idx < end; ++idx) {
                std::size_t row = cData.row_index[idx]; // Row index
                T value = cData.values[idx];           // Non-zero value

                // Insert into COOmap (uData)
                uData[{row, col}] = value;
            }
        }
    }
    
    cData.values.clear();
    cData.row_index.clear();
    cData.col_index.clear();
    compressed = false; // Mark as uncompressed
}

// METHOD TO PRINT THE MATRIX (WORKS BOTH FOR COMPRESSED AND UNCOMPRESSED STATE)
template <typename T, StorageOrder Order>
void Matrix<T, Order>::print() const{
    if constexpr (Order == ROW_MAJOR){
        // we make a distinction between arithmetic types and complex numbers
        if constexpr (std::is_same_v<T, std::complex<double>>) {
            if(compressed){
                std::cout << "Matrix is compressed" << std::endl;
                std::cout << "Matrix dimensions : (" << rows << ", " << cols << ")" << std::endl;
                std::cout << "Number of non-zero values : " << nnz << std::endl;
                // Row index
                std::cout << "Row index : " << std::endl;
                for(std::size_t i = 0; i < cData.row_index.size(); ++i){
                    std::cout << cData.row_index[i] << " ";
                }
                std::cout << std::endl;
                // Column index
                std::cout << "Column index : " << std::endl;
                for(std::size_t i = 0; i < cData.col_index.size(); ++i){
                    std::cout << cData.col_index[i] << " ";
                }
                std::cout << std::endl;
                // Values
                std::cout << "Values : " << std::endl;
                for(std::size_t i = 0; i < cData.values.size(); ++i){
                    std::cout << cData.values[i].real() << "+i" << cData.values[i].imag() << " ";
                }
                std::cout << std::endl;
            } else{
                std::cout << "Matrix is uncompressed" << std::endl;
                std::cout << "Matrix dimensions : (" << rows << ", " << cols << ")" << std::endl;
                std::cout << "Number of non-zero values : " << nnz << std::endl;
                // Uncompressed data
                for(auto& [key, value] : uData){
                    std::cout << "[" << key.first << ", " << key.second << "] : " << std::fixed << std::setprecision(2) <<
                    value.real() << "+i" << value.imag() << std::endl;
                }
            }
        } else {
            if(compressed){
                std::cout << "Matrix is compressed" << std::endl;
                std::cout << "Matrix dimensions : (" << rows << ", " << cols << ")" << std::endl;
                std::cout << "Number of non-zero values : " << nnz << std::endl;
                // Row index
                std::cout << "Row index : " << std::endl;
                for(std::size_t i = 0; i < cData.row_index.size(); ++i){
                    std::cout << cData.row_index[i] << " ";
                }
                std::cout << std::endl;
                // Column index
                std::cout << "Column index : " << std::endl;
                for(std::size_t i = 0; i < cData.col_index.size(); ++i){
                    std::cout << cData.col_index[i] << " ";
                }
                std::cout << std::endl;
                // Values
                std::cout << "Values : " << std::endl;
                for(std::size_t i = 0; i < cData.values.size(); ++i){
                    std::cout << cData.values[i] << " ";
                }
                std::cout << std::endl;
            } else{
                std::cout << "Matrix is uncompressed" << std::endl;
                std::cout << "Matrix dimensions : (" << rows << ", " << cols << ")" << std::endl;
                std::cout << "Number of non-zero values : " << nnz << std::endl;
                // Uncompressed data
                for(auto& [key, value] : uData){
                    std::cout << "[" << key.first << ", " << key.second << "] : " << std::fixed << std::setprecision(2) << value << std::endl;
                }
            }
        }       
    } else {
        // we make a distinction between arithmetic types and complex numbers
        if constexpr (std::is_same_v<T, std::complex<double>>)
            if(compressed){
                std::cout << "Matrix is compressed" << std::endl;
                std::cout << "Matrix dimensions : (" << cols << ", " << rows << ")" << std::endl;
                std::cout << "Number of non-zero values : " << nnz << std::endl;
                // Column index
                std::cout << "Column index : " << std::endl;
                for(std::size_t i = 0; i < cData.col_index.size(); ++i){
                    std::cout << cData.col_index[i] << " ";
                }
                std::cout << std::endl;
                // Row index
                std::cout << "Row index : " << std::endl;
                for(std::size_t i = 0; i < cData.row_index.size(); ++i){
                    std::cout << cData.row_index[i] << " ";
                }
                std::cout << std::endl;
                // Values
                std::cout << "Values : " << std::endl;
                for(std::size_t i = 0; i < cData.values.size(); ++i){
                    std::cout << cData.values[i].real() << "+i" << cData.values[i].imag() << " ";
                }
                std::cout << std::endl;
            } else{
                std::cout << "Matrix is uncompressed" << std::endl;
                std::cout << "Matrix dimensions : (" << rows << ", " << cols << ")" << std::endl;
                std::cout << "Number of non-zero values : " << nnz << std::endl;
                // Uncompressed data
                for(auto& [key, value] : uData){
                    std::cout << "[" << key.first << ", " << key.second << "] : " << std::fixed << std::setprecision(2) <<
                    value.real() << "+i" << value.imag() << std::endl;
                }
            }
        else {
            if(compressed){
                std::cout << "Matrix is compressed" << std::endl;
                std::cout << "Matrix dimensions : (" << cols << ", " << rows << ")" << std::endl;
                std::cout << "Number of non-zero values : " << nnz << std::endl;
                // Column index
                std::cout << "Column index : " << std::endl;
                for(std::size_t i = 0; i < cData.col_index.size(); ++i){
                    std::cout << cData.col_index[i] << " ";
                }
                std::cout << std::endl;
                // Row index
                std::cout << "Row index : " << std::endl;
                for(std::size_t i = 0; i < cData.row_index.size(); ++i){
                    std::cout << cData.row_index[i] << " ";
                }
                std::cout << std::endl;
                // Values
                std::cout << "Values : " << std::endl;
                for(std::size_t i = 0; i < cData.values.size(); ++i){
                    std::cout << cData.values[i] << " ";
                }
                std::cout << std::endl;
            } else{
                std::cout << "Matrix is uncompressed" << std::endl;
                std::cout << "Matrix dimensions : (" << rows << ", " << cols << ")" << std::endl;
                std::cout << "Number of non-zero values : " << nnz << std::endl;
                // Uncompressed data
                for(auto& [key, value] : uData){
                    std::cout << "[" << key.first << ", " << key.second << "] : " << std::fixed << std::setprecision(2) << value << std::endl;
                }
            }
        }
    }
}

// METHOD TO PRINT THE MATRIX IN DENSE FORMAT
template <typename T, StorageOrder Order>
void Matrix<T, Order>::printdense() const{

    // Handle the case when the matrix is too big to print in dense format (50 is arbitrary)
    if(rows > MAX_SIZE_PRINTABLE || cols > MAX_SIZE_PRINTABLE){
        std::cerr << "Matrix is too big to print in dense format" << std::endl;
        return;
    }

    std::cout << "Matrix in dense format:" << std::endl;
    for(std::size_t i = 0; i < rows; ++i) {
        for(std::size_t j = 0; j < cols; ++j) {
            // We use the operator() to access the values
            std::cout << this->operator()(i, j) << " ";
        }
        std::cout << std::endl;
    }

}

// METHOD TO COMPUTE THE NORM OF THE MATRIX
template<typename T, StorageOrder Order>
template<NormType N>
double Matrix<T, Order>::norm() const{
    // First we make sure that the matrix contains coefficients that are either arithmetic or std::complex
    // This is fundamental, otherwise we can't compute the norm
    static_assert(std::is_arithmetic_v<T> || std::is_same_v<T, std::complex<double>>,
        "Matrix norm is only defined for arithmetic or std::complex types");

    // We define a lambda function to compute the absolute value of the coefficients
    // This works for both arithmetic and std::complex types
    auto abs_val = [](const T& val) -> double{
        return std::abs(val);
    };

    // Now we can start computing the norms
    // We make a distinction between the Norm Types first
    
           
    if constexpr ( N == FROBENIUS) {
        // This is independent of the StorageOrder
        // We simply iterate over the data structures and sum the squares of all the values,
        // then return the squared root
        if (compressed) {
            double result = 0;
            for (auto & value : cData.values){
                result += abs_val(value)*abs_val(value);
            }
            return std::sqrt(result); 
        } else {
            double result = 0;
            for (auto & [key, value] : uData) {
                result += abs_val(value)*abs_val(value);
            }
            return std::sqrt(result);
        }
    } else if constexpr ( N == ONE ) {
        if constexpr (Order == ROW_MAJOR) {
            if (compressed) {
                // We create an auxiliary vector to store the sums on each column
                std::vector<double> column_sums(cols, 0.0);
                // We now loop over all values and sum their absolute values to the corresponding slot
                // in the auxiliary vector
                for (std::size_t i = 0; i < nnz; ++i) {
                    column_sums[cData.col_index[i]] += abs_val(cData.values[i]);
                }
                // We return the maximum value among the ones in the auxiliary vector
                return *std::max_element(column_sums.begin(), column_sums.end());
            } else { // UNCOMPRESSED
                // We create an auxiliary vector to store the sums on each column
                std::vector<double> column_sums(cols, 0.0);
                // We now loop over all values and sum their absolute values to the corresponding slot
                // in the auxiliary vector
                for (auto & [key, value] : uData) {
                    column_sums[key.second] += abs_val(value);
                }
                // We return the maximum value among the ones in the auxiliary vector
                return *std::max_element(column_sums.begin(), column_sums.end());
            }
        } else { // COL_MAJOR
            if (compressed){
                // We iterate over the col_index vector to get the number of elements from the current column
                // and then extract their values. We perform a comparison with the current column when changing it
                // to see if we have a new maximum
                double result = 0.0;
                double sum;
                std::size_t elementsPerColumn;
                std::size_t firstElem = 0;
                for (std::size_t i = 1; i < cols + 1; ++i) {
                    elementsPerColumn = cData.col_index[i] - cData.col_index[i-1];
                    sum = 0;
                    for (std::size_t j = 0; j < elementsPerColumn; ++j) {
                        sum += abs_val(cData.values[firstElem + j]);
                    }
                    firstElem += elementsPerColumn;
                    result = std::max(result, sum);
                }
                return result;
            } else { // UNCOMPRESSED
                // In this case we iterate through the ordered map exploiting the fact that for a COL_MAJOR matrix
                // this is automatically done by columns. We perform a check when we change column to see if we have a
                // new maximum
                double result = 0.0;
                double sum = 0;
                std::size_t currColumn = uData.begin()->first.second;
                for (auto & [key, value] : uData) {
                    if (key.second == currColumn) {
                        sum += abs_val(value); 
                    } else {
                        result = std::max(result, sum);
                        currColumn = key.second;
                        sum = abs_val(value);
                    }
                }
                // we have to perform a check here too in the case the maximum was in the last column
                result = std::max(result, sum);
                return result;
            }
        }
    } else if constexpr ( N == INFTY) {
        if constexpr ( Order == ROW_MAJOR) {
            if (compressed){
                // We iterate over the row_index vector to get the number of elements from the current row
                // and then extract their values. We perform a comparison with the current row when changing it
                // to see if we have a new maximum
                double result = 0.0;
                double sum;
                std::size_t elementsPerRow;
                std::size_t firstElem = 0;
                for (std::size_t i = 1; i < rows + 1; ++i) {
                    elementsPerRow = cData.row_index[i] - cData.row_index[i-1];
                    sum = 0;
                    for (std::size_t j = 0; j < elementsPerRow; ++j) {
                        sum += abs_val(cData.values[firstElem + j]);
                    }
                    firstElem += elementsPerRow;
                    result = std::max(result, sum);
                }
                return result;
            } else {
                // In this case we iterate through the ordered map exploiting the fact that for a ROW_MAJOR matrix
                // this is automatically done by rows. We perform a check when we change row to see if we have a
                // new maximum
                double result = 0.0;
                double sum = 0;
                std::size_t currRow = uData.begin()->first.first;
                for (auto & [key, value] : uData) {
                    if (key.first == currRow) {
                        sum += abs_val(value); 
                    } else {
                        result = std::max(result, sum);
                        currRow = key.first;
                        sum = abs_val(value);
                    }
                }
                // we have to perform a check here too in the case the maximum was in the last row
                result = std::max(result, sum);
                return result;
            }
        } else { // COL_MAJOR
            if (compressed) {
                // We create an auxiliary vector to store the sums on each row
                std::vector<double> row_sums(rows, 0.0);
                // We now loop over all values and sum their absolute values to the corresponding slot
                // in the auxiliary vector
                for (std::size_t i = 0; i < nnz; ++i) {
                    row_sums[cData.row_index[i]] += abs_val(cData.values[i]);
                }
                // We return the maximum value among the ones in the auxiliary vector
                return *std::max_element(row_sums.begin(), row_sums.end());
            } else {    // UNCOMPRESSED
                // We create an auxiliary vector to store the sums on each row
                std::vector<double> row_sums(cols, 0.0);
                // We now loop over all values and sum their absolute values to the corresponding slot
                // in the auxiliary vector
                for (auto & [key, value] : uData) {
                    row_sums[key.first] += abs_val(value);
                }
                // We return the maximum value among the ones in the auxiliary vector
                return *std::max_element(row_sums.begin(), row_sums.end());
            }
        }
    }
}

// OVERLOAD OPERATOR * FOR MATRIX-VECTOR MULTIPLICATION
template<typename T, StorageOrder Order>
std::vector<T> algebra::operator*(const Matrix<T, Order>& M, const std::vector<T>& v){
    // Check if the vector is compatible with the matrix
    if(M.getCols() != v.size()){
        std::cerr << "Error: Matrix and vector dimensions do not match" << std::endl;
        return {};
    }

    // Create the result vector to fill
    std::vector<T> result(M.getRows(), 0);

    if(M.is_compressed()){

        if constexpr(Order == ROW_MAJOR){
            // Iterate over the rows of the compressed matrix
            for (std::size_t row = 0; row < M.getRows(); ++row) {
                std::size_t start = M.cData.row_index[row];
                std::size_t end = M.cData.row_index[row + 1];

                // Iterate over the non-zero elements in the current row
                for (std::size_t idx = start; idx < end; ++idx) {
                    std::size_t col = M.cData.col_index[idx]; // Column index
                    T value = M.cData.values[idx];           // Non-zero value

                    // Multiply and accumulate
                    result[row] += value * v[col];
                }
            }    
        }
        else{
            // Iterate over the columns of the compressed matrix
            for (std::size_t col = 0; col < M.getCols(); ++col) {
                std::size_t start = M.cData.col_index[col];
                std::size_t end = M.cData.col_index[col + 1];

                // Iterate over the non-zero elements in the current column
                for (std::size_t idx = start; idx < end; ++idx) {
                    std::size_t row = M.cData.row_index[idx]; // Row index
                    T value = M.cData.values[idx];           // Non-zero value

                    // Multiply and accumulate
                    result[row] += value * v[col];
                }
            }
        }

    }else{
         // Iterate over the uncompressed data
         for (const auto& [key, value] : M.uData) {
            result[key.first] += value * v[key.second];
        }
    }

    // Return the result vector
    return result;
}

// OVERLOAD OPERATOR * FOR MATRIX-MATRIX MULTIPLICATION
template<typename T, StorageOrder Order>
Matrix<T,Order> algebra::operator*(const Matrix<T, Order>& M1, const Matrix<T, Order>& M2){
    // Check if the matrices are compatible for multiplication
    if(M1.getCols() != M2.getRows()){
        std::cerr << "Error: Matrix dimensions do not match for multiplication" << std::endl;
        return {};
    }
    
    // Create a result matrix
    Matrix<T,Order> result;
    result.resize(M1.getRows(), M2.getCols());

    if constexpr(Order == ROW_MAJOR){
        // In order to make the multiplication easier, we transpose the second matrix (we create a copy of it so the original one is not modified)
        Matrix<T,Order> M2_temp = M2;
        Matrix<T,Order> M1_temp = M1;
        M2_temp.transpose();

        // If at least one of the matrices is compressed we compress the other one and the result will be a compressed matrix
        if(M1.is_compressed() || M2.is_compressed()){
            M1_temp.compress();
            M2_temp.compress();

            // We iterate over the rows of the first matrix
            std::size_t rows1 = M1_temp.getRows();
            std::size_t rows2 = M2_temp.getRows();
            for(std::size_t i=0; i<rows1; ++i){
                // We compute the number of elements in the current row of the first matrix
                std::size_t row1 = M1_temp.cData.row_index[i+1] - M1_temp.cData.row_index[i];
                // We iterate over the rows of the second matrix
                for(std::size_t j=0; j<rows2; ++j){
                    // We compute the number of elements in the current row of the second matrix
                    std::size_t row2 = M2_temp.cData.row_index[j+1] - M2_temp.cData.row_index[j];
                    std::size_t counter1 = 0;
                    std::size_t counter2 = 0;
                    T sum = 0;
                    // We iterate over the elements of both rows, storing the result only if the columns index is the same
                    while(counter1<row1 && counter2<row2){
                        if(M1_temp.cData.col_index[M1_temp.cData.row_index[i]+counter1] == M2_temp.cData.col_index[M2_temp.cData.row_index[j]+counter2]){
                            sum += M1_temp.cData.values[M1_temp.cData.row_index[i]+counter1] * M2_temp.cData.values[M2_temp.cData.row_index[j]+counter2];
                            ++counter1;
                            ++counter2;
                        }
                        else if(M1_temp.cData.col_index[M1_temp.cData.row_index[i]+counter1] < M2_temp.cData.col_index[M2_temp.cData.row_index[j]+counter2]){
                            ++counter1;
                        }
                        else if(M1_temp.cData.col_index[M1_temp.cData.row_index[i]+counter1] > M2_temp.cData.col_index[M2_temp.cData.row_index[j]+counter2]){
                            ++counter2;
                        }
                    }
                    if(sum != 0){
                        result.set(i,j,sum);
                    }
                }
            }
            // We compress the result
            result.compress();
        }
        // If both matrices are uncompressed, then the result will be an uncompressed matrix
        else{
            // We iterate over all the rows of the first matrix
            std::size_t rows1 = M1_temp.getRows();
            std::size_t rows2 = M2_temp.getRows();
            auto bk1 = M1_temp.uData.cbegin();
            for(std::size_t i=0; i<rows1; ++i){
                auto it2 = M2_temp.uData.cbegin();
                // We iterate over all the rows of the second matrix
                for(std::size_t j=0; j<rows2; ++j){
                    auto it1 = bk1;
                    T sum = 0;
                    // If we find the same column index in both rows, we compute the multiplication and store it in the result matrix
                    while(it1!=M1_temp.uData.cend() && it2!=M2_temp.uData.cend() && it1->first.first==i && it2->first.first==j){
                        if(it1->first.second == it2->first.second){
                            sum += it1->second * it2->second;
                            ++it1;
                            ++it2;
                        }
                        else if(it1->first.second < it2->first.second){
                            ++it1;
                        }
                        else if(it1->first.second > it2->first.second){
                            ++it2;
                        }
                    }
                    // After the while loop, we make sure that the next element is not in the same row of the second matrix
                    while(it2!=M2_temp.uData.cend() && it2->first.first==j){
                        ++it2;
                    }
                    if(sum != 0){
                        result.set(i,j,sum);
                    }
                }
                // We make sure that the next element is not in the same row of the first matrix
                while(bk1!=M1_temp.uData.cend() && bk1->first.first==i){
                    ++bk1;
                }
            }
        }

        return result;
    }
    if constexpr(Order == COL_MAJOR){
        // In order to make the multiplication easier, we transpose the first matrix (we create a copy of it so the original one is not modified)
        Matrix<T,Order> M2_temp = M2;
        Matrix<T,Order> M1_temp = M1;
        M1_temp.transpose();

        // If at least one of the matrices is compressed we compress the other one and the result will be a compressed matrix
        if(M1.is_compressed() || M2.is_compressed()){
            M1_temp.compress();
            M2_temp.compress();

            // We iterate over the rows of the first matrix (cols of the transposed matrix)
            std::size_t cols1 = M1_temp.getCols();
            std::size_t cols2 = M2_temp.getCols();
            for(std::size_t i=0; i<cols1; ++i){
                // We compute the number of elements in the current column of the first (transposed) matrix
                std::size_t col1 = M1_temp.cData.col_index[i+1] - M1_temp.cData.col_index[i];
                // We iterate over the columns of the second matrix
                for(std::size_t j=0; j<cols2; ++j){
                    // We compute the number of elements in the current column of the second matrix
                    std::size_t col2 = M2_temp.cData.col_index[j+1] - M2_temp.cData.col_index[j];
                    std::size_t counter1 = 0;
                    std::size_t counter2 = 0;
                    T sum = 0;
                    // We iterate over the elements of both columns, storing the result only if the row index is the same
                    while(counter1<col1 && counter2<col2){
                        if(M1_temp.cData.row_index[M1_temp.cData.col_index[i]+counter1] == M2_temp.cData.row_index[M2_temp.cData.col_index[j]+counter2]){
                            sum += M1_temp.cData.values[M1_temp.cData.col_index[i]+counter1] * M2_temp.cData.values[M2_temp.cData.col_index[j]+counter2];
                            ++counter1;
                            ++counter2;
                        }
                        else if(M1_temp.cData.row_index[M1_temp.cData.col_index[i]+counter1] < M2_temp.cData.row_index[M2_temp.cData.col_index[j]+counter2]){
                            ++counter1;
                        }
                        else if(M1_temp.cData.row_index[M1_temp.cData.col_index[i]+counter1] > M2_temp.cData.row_index[M2_temp.cData.col_index[j]+counter2]){
                            ++counter2;
                        }
                    }
                    if(sum != 0){
                        result.set(i,j,sum);
                    }
                }
            }
            // We compress the result
            result.compress();
        }
        // If both matrices are uncompressed, then the result will be an uncompressed matrix
        else{
            // We iterate over all the columns of the first (transposed) matrix
            std::size_t cols1 = M1_temp.getCols();
            std::size_t cols2 = M2_temp.getCols();
            auto bk1 = M1_temp.uData.cbegin();
            for(std::size_t i=0; i<cols1; ++i){
                auto it2 = M2_temp.uData.cbegin();
                // We iterate over all the columns of the second matrix
                for(std::size_t j=0; j<cols2; ++j){
                    auto it1 = bk1;
                    T sum = 0;
                    // If we find the same row index in both columns, we compute the multiplication and store it in the result matrix
                    while(it1!=M1_temp.uData.cend() && it2!=M2_temp.uData.cend() && it1->first.second==i && it2->first.second==j){
                        if(it1->first.first == it2->first.first){
                            sum += it1->second * it2->second;
                            ++it1;
                            ++it2;
                        }
                        else if(it1->first.first < it2->first.first){
                            ++it1;
                        }
                        else if(it1->first.first > it2->first.first){
                            ++it2;
                        }
                    }
                    // After the while loop, we make sure that the next element is not in the same column of the second matrix
                    while(it2!=M2_temp.uData.cend() && it2->first.second==j){
                        ++it2;
                    }
                    if(sum != 0){
                        result.set(i,j,sum);
                    }
                }
                // We make sure that the next element is not in the same column of the first matrix
                while(bk1!=M1_temp.uData.cend() && bk1->first.second==i){
                    ++bk1;
                }
            }
        }

        return result;
    }
}