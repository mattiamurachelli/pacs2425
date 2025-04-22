#include "Minimizer.hpp"
#include <iostream>

Minimizer::Minimizer(const std::function<double(vectorExp)> &f,
             	     const std::function<vectorExp(vectorExp)> &grad,
                     vectorExp x0,
                     double ts, double tr, double alpha_0,   
                     std::size_t mi)
: f(f),grad(grad),x(x0),tol_step(ts),tol_res(tr),alpha_k(alpha_0),max_it(mi),it(0), alpha_0(alpha_0){};

const vectorExp& Minimizer::get_x() const{
	return x;
}

const std::size_t& Minimizer::get_iter() const{
	return it;
}

bool Minimizer::converged(double step, double residual) const {
	bool flag = false;

	if(std::abs(step)<tol_step || std::abs(residual)<tol_res || it >= max_it)
		flag = true;

	return flag;
}

void Minimizer::solve(){

	double res  = tol_res + 1;
	double step = tol_step + 1;

	vectorExp xk = x;
	vectorExp gradient({0,0});

	while(!converged(step,res)){

		// Compute alpha
		compute_alpha();

		// Compute gradient of f in xk
		gradient = grad(xk);

		// gradient clipping for very steep functions
		if (gradient.norm() > 1){
			gradient = gradient*(1/gradient.norm());
		}

		x = xk - gradient*alpha_k;

		// update res and step
		res = grad(x).norm();
		step = (x-xk).norm();

		// update xk
		xk = x;

		// Update iterations number
		++it;
		
	}

}
