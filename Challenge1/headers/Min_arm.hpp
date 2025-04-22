#ifndef MIN_ARM_HPP
#define MIN_ARM_HPP

#include "vectorExp.hpp"
#include "Minimizer.hpp"
#include <functional>

class Min_arm : public Minimizer {
    private :
        double sigma;
        void compute_alpha() override;

    public :
        Min_arm(const std::function<double(vectorExp)> &f,
                const std::function<vectorExp(vectorExp)> &grad,
                double sigma,vectorExp x0=vectorExp({0,0}),
                double ts=1e-6, double tr=1e-6, double alpha_0=1,
                std::size_t mi=100);
};

#endif // MIN_ARM_HPP
