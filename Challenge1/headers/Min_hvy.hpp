#ifndef MIN_HVY_HPP
#define MIN_HVY_HPP

#include "vectorExp.hpp"
#include "Minimizer.hpp"
#include <functional>

class Min_hvy : public Minimizer {
    private :
        double eta;
        double mu;
        vectorExp delta{0,0};
        vectorExp xk{0,0};
        void compute_alpha() override;

    public :
        Min_hvy(const std::function<double(vectorExp)> &f,
                const std::function<vectorExp(vectorExp)> &grad,
                double eta, double mu, vectorExp x0=vectorExp({0,0}),
                double ts=1e-6, double tr=1e-6, double alpha_0=1,
                std::size_t mi=100);

        void solve() override;
};

#endif // MIN_EXP_HPP