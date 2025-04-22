#ifndef MIN_INV_HPP
#define MIN_INV_HPP

#include "vectorExp.hpp"
#include "Minimizer.hpp"
#include <functional>

class Min_inv : public Minimizer {
    private :
        double mu;
        void compute_alpha() override;

    public:
        Min_inv(const std::function<double(vectorExp)> &f,
                const std::function<vectorExp(vectorExp)> &grad,
                double mu, vectorExp x0=vectorExp({0,0}),
                double ts=1e-6, double tr=1e-6, double alpha_0=1,
                std::size_t mi=100);
};

#endif //MIN_INV_HPP