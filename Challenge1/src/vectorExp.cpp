#include "vectorExp.hpp"
#include <cmath>
#include <iostream>
#include <limits>

double vectorExp::norm() const {
    if(empty()){
        std::cerr << "Norm undefined: empty vector" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum = 0;
    for(auto e : data){
        sum += e*e;
    }
    return std::sqrt(sum);
}

vectorExp vectorExp::operator*(const double c) const {
    vectorExp res{};

    if(empty()){
        std::cerr << "Undefined result for *: empty vector" << std::endl;
        return res;
    }

    for(auto e : data)
        res.data.emplace_back(c*e);

    return res;
}

double vectorExp::operator[](const std::size_t i) const {
    if(empty()){
        std::cerr << "Undefined result for []: empty vector" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    if(i >= size()){
        std::cerr << "Index out of range" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    return data[i];
}

std::size_t vectorExp::size() const {
    return data.size();
}

void vectorExp::print() const {
    if(empty()){
        std::cout << "Empty vector" << std::endl;
        return;
    }

    std::cout << "[ ";
    for(auto e : data){
        std::cout << e << " ";
    }
    std::cout << "]" << std::endl;

    return;
}

bool vectorExp::empty() const {
    return data.empty();
}

vectorExp operator+(const vectorExp & v, const vectorExp & w) {
    vectorExp res{};
    
    if(v.empty() || w.empty()){
        std::cerr << "Undefined result for +: one operand is empty" << std::endl;
        return res;
    }

    if(v.size()!=w.size()) {
        std::cerr << "Invalid + operation: size mismatch" << std::endl;
        return res;
    }

    for(std::size_t i=0; i<w.size(); i++)
        res.data.emplace_back(v[i]+w[i]);
    
    return res;
}

vectorExp operator-(const vectorExp & v, const vectorExp & w) {
    vectorExp res{};

    if(v.empty() || w.empty()){
        std::cerr << "Undefined result for -: one operand is empty" << std::endl;
        return res;
    }

    if(v.size()!=w.size()) {
        std::cerr << "Invalid - operation: size mismatch" << std::endl;
        return res;
    }

    for(std::size_t i=0; i<w.size(); i++)
        res.data.emplace_back(v[i]-w[i]);
    
    return res;
}

vectorExp operator*(const double c, const vectorExp & v) {
    vectorExp res{};

    if(v.empty()){
        std::cerr << "Undefined result for *: empty vector" << std::endl;
        return res;
    }

    for(auto e : v.data)
        res.data.emplace_back(c*e);

    return res;   
}

double vectorExp::dot(const vectorExp &v) const {
    if(empty() || v.empty()){
        std::cerr << "Undefined result for dot: one operand is empty" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    if(size()!=v.size()){
        std::cerr << "Invalid dot operation: size mismatch" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum = 0;
    for(std::size_t i=0; i<size(); i++)
        sum += data[i]*v[i];

    return sum; 
}