#include "Min_arm.hpp"
#include <iostream>

Min_arm::Min_arm(const std::function<double(vectorExp)> &f,
                const std::function<vectorExp(vectorExp)> &grad,
                double sigma, vectorExp x0,
                double ts, double tr, double alpha_0,
                std::size_t mi) :
                Minimizer(f,grad,x0,ts,tr,alpha_0,mi), sigma(sigma) {};

void Min_arm::compute_alpha(){

	double alpha_temp=alpha_0;
	int iterations = 0;
	vectorExp gradient = grad(x);
	double gradient_norm = gradient.norm();
	double f_x = f(x);
	double lhs = f_x - f(x-gradient*alpha_temp);
	double rhs = sigma*alpha_temp*gradient_norm*gradient_norm;

	while(lhs<rhs && iterations<25){
		++iterations;
		alpha_temp /= 2;
		lhs = f_x - f(x-gradient*alpha_temp);
		rhs = sigma*alpha_temp*gradient_norm*gradient_norm;

		std::cout << "Value for alpha_temp " << alpha_temp << std::endl;

	}

	alpha_k = alpha_temp;
}
