#ifndef HEADER_H_MINIMIZER
#define HEADER_H_MINIMIZER

#include "vectorExp.hpp"
#include <functional>

class Minimizer{
protected:
	const std::function<double(vectorExp)> &f;    	  // Function to minimize
	const std::function<vectorExp(vectorExp)> &grad;  // Gradient of the function
	vectorExp x;							  // Solution
	double tol_step;						  // Tolerance on the step
	double tol_res;							  // Tolerance on the residual
	double alpha_k;							  // Step alpha
	std::size_t max_it;					      // Maximum number of iterations
	std::size_t it;					          // Number of iterations
	double alpha_0;							  // Initial value for alpha

	// Check if converged
	bool converged(double step, double residual) const;

	// Compute alpha (virtual)
	virtual void compute_alpha() = 0;
	
public:
	Minimizer(const std::function<double(vectorExp)> &f,
			  const std::function<vectorExp(vectorExp)> &grad, 
			  vectorExp x0,
			  double ts, double tr, double alpha_0,
			  std::size_t mi);

			  virtual ~Minimizer() = default;
			  
	// Getters
	const vectorExp&    get_x()    const; 
	const std::size_t&  get_iter() const;

	// Solver function
	virtual void solve();
	
};

#endif
