#include "Min_hvy.hpp"
#include <cmath>
#include <iostream>

Min_hvy::Min_hvy(const std::function<double(vectorExp)> &f,
                const std::function<vectorExp(vectorExp)> &grad,
                double eta, double mu, vectorExp x0,
                double ts, double tr, double alpha_0,
                std::size_t mi) : Minimizer(f, grad, x0, ts, tr, alpha_0, mi), eta(eta), mu(mu) {}

void Min_hvy::compute_alpha() {
    // we were not able to find a suitable compute_alpha() for this method
    // as the function is too steep near the minimum
    
    // alpha_k = alpha_0; (constant step size)
    // alpha_k = alpha_0*std::exp(-mu*it); (exponential decay)
    // alpha_k = alpha_0 / (1 + mu*it); (inverse decay)

    // Barzilai-Borwein method
    if(it == 0){
        alpha_k = alpha_0;
    } else {
        vectorExp s = x - xk;
        vectorExp y = grad(x) - grad(xk);
        alpha_k = s.dot(s) / s.dot(y);

        double numerator = s.dot(s);
        double denominator = s.dot(y);

        // make sure we are not dividing by zero
        if (std::abs(denominator) > 1e-8) {
            alpha_k = numerator / denominator;
        } else {
            alpha_k = alpha_0;
        }

        // Clamp the value of alpha_k
        alpha_k = std::min(std::max(alpha_k, 1e-4), 1e-2);
    }
}

void Min_hvy::solve() {

    double res  = tol_res + 1;
	double step = tol_step + 1;

	vectorExp xk = x;
	vectorExp gradient = grad(xk);
    delta = -alpha_0*gradient;

	while(!converged(step,res)){

		// Compute alpha
		compute_alpha();
        
        // update pass
		x = xk + delta;
        gradient = grad(x);
        delta = eta*delta - alpha_k*gradient;

		// update res and step
		res = gradient.norm();
		step = (x-xk).norm();

		// update xk
		xk = x;

		// Update iterations number
		++it;

        // output print
        std::cout << "Current value of x : ";
        x.print();
        std::cout << "Current value of f(x) : " << f(x) << std::endl;
	}
}
