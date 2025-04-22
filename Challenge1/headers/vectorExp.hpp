#ifndef VECTOREXP_HPP
#define VECTOREXP_HPP

#include <vector>
#include <initializer_list>

class vectorExp {
    friend vectorExp operator+(const vectorExp &, const vectorExp &);
    friend vectorExp operator-(const vectorExp &, const vectorExp &); 
    friend vectorExp operator*(const double, const vectorExp &);

    std::vector<double> data{};

    public:
    vectorExp() = default;
    vectorExp(std::initializer_list<double> values): data(values){};
    double norm() const;
    vectorExp operator*(const double c) const;
    double operator[](const std::size_t i) const;
    std::size_t size() const;
    void print() const;
    bool empty() const;
    double dot(const vectorExp & v) const;

};

vectorExp operator+(const vectorExp & v, const vectorExp & w);
vectorExp operator-(const vectorExp & v, const vectorExp & w);
vectorExp operator*(const double, const vectorExp &);

#endif //VECTOREXP_HPP