#include "Min_inv.hpp"
#include <cmath>

Min_inv::Min_inv(const std::function<double(vectorExp)> &f,
                const std::function<vectorExp(vectorExp)> &grad,
                double mu, vectorExp x0,
                double ts, double tr, double alpha_0,
                std::size_t mi) : Minimizer(f, grad, x0, ts, tr, alpha_0, mi), mu(mu) {}
                
                
void Min_inv::compute_alpha() {
    alpha_k = alpha_0 / (1 + mu*it);
}