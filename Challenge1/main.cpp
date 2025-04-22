// Main file of our solution for the 1st challenge of the PACS Course

#include <iostream>
#include <functional>
#include <set>
#include "vectorExp.hpp"
#include "Minimizer.hpp"
#include "Min_exp.hpp"
#include "Min_inv.hpp"
#include "Min_arm.hpp"
#include "Min_hvy.hpp"
#include <limits>

using Fun = std::function<double(vectorExp)>;
using Grad = std::function<vectorExp(vectorExp)>; 

double function(vectorExp x) {
    return 2*x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + 2*x[1]*x[1] + 2*x[0];
    // return (x[0]-1)*(x[0]-1) + (x[1]-1)*(x[1]-1);
}

vectorExp gradient(vectorExp x) {
    return vectorExp({2*x[1] + 16*x[0]*x[0]*x[0] + 2, 2*x[0] + 4*x[1]});
    // return vectorExp({2*(x[0]-1), 2*(x[1]-1)});
}

int main() {
    
    std::set methods = {0, 1, 2, 3};
    int user_choice = -1;
    std::string custom_settings = "N";

    Fun f = function;
    Grad gradf = gradient;

    // declare a pointer to the base class in order to exploit polymorphism
    Minimizer * minimizer = nullptr;

    double mu = 0.2;
    double sigma = 0.25;
    double eta = 0.9;

    // list here the implemented methods
    std::cout << "[0] Exponential decay rule : " << std::endl;
    std::cout << "[1] Inverse decay rule : " << std::endl;
    std::cout << "[2] Armijo rule : " << std::endl;
    std::cout << "[3] Heavy-ball rule : " << std::endl;
    // ask the user the method he wants to use
    while(methods.find(user_choice) == methods.end()){
        std::cout << "Please, choose the optimization method you want to use : ";
        std::cin >> user_choice;
    }
    // and also ask him if he wants to use custom settings or not
    std::cout << "Do you want to use custom settings for the optimization method? (Y/N) : ";
    std::cin >> custom_settings;
    if(custom_settings == "Y" || custom_settings == "y" || custom_settings == "Yes" || custom_settings == "yes"){

        double x_coord = 0, y_coord = 0;
        std::cout << "Insert x coordinate for the starting point : ";
        std::cin >> x_coord;
        std::cout << "Insert y coordinate for the starting point : ";
        std::cin >> y_coord;
        vectorExp x0({x_coord, y_coord});

        double tolerance = std::numeric_limits<double>::epsilon()*1000;
        std::cout << "Insert value for tolerance : ";
        std::cin >> tolerance;

        double alpha_initial = 1;
        std::cout << "Insert initial value for alpha : ";
        std::cin >> alpha_initial;

        std::size_t max_it = 100;
        std::cout << "Insert maximum number of iterations : ";
        std::cin >> max_it;

         // instantiate the solver using the correct one based on the user's choice and the custom settings 
        if(user_choice == 0){
            std::cout << "Please insert the value for mu (suggested value = 0.2) : ";
            std::cin >> mu;
            minimizer = new Min_exp(f, gradf, mu, x0, tolerance, tolerance, alpha_initial, max_it);
        } else if(user_choice == 1){
            std::cout << "Please insert the value for mu (suggested value = 0.2) : ";
            std::cin >> mu;
            minimizer = new Min_inv(f, gradf, mu, x0, tolerance, tolerance, alpha_initial, max_it);
        } else if(user_choice == 2){
            std::cout << "Please insert the value for sigma (suggested value = 0.25) : ";
            std::cin >> sigma;
            minimizer = new Min_arm(f, gradf, sigma, x0, tolerance, tolerance, alpha_initial, max_it);
        } else if(user_choice == 3) {
            std::cout << "Please insert the value for eta (suggested value = 0.9) : ";
            std::cin >> eta;
            std::cout << "Please insert the value for mu (suggested value = 0.2) : ";
            std::cin >> mu;
            minimizer = new Min_hvy(f, gradf, eta, mu, x0, tolerance, tolerance, alpha_initial, max_it);
        }
    } else {
        if(custom_settings != "N" && custom_settings != "n" && custom_settings != "No" && custom_settings != "no"){
            std::cout << "Failed to provide a valid input, default settings will be used." << std::endl;
        }
        // instantiate the solver using the correct one based on the user's choice and default settings
        if(user_choice == 0){
            std::cout << "Please insert the value for mu (suggested value = 0.2) : ";
            std::cin >> mu;
            minimizer = new Min_exp(f, gradf, mu);
        } else if(user_choice == 1){
            std::cout << "Please insert the value for mu (suggested value = 0.2) : ";
            std::cin >> mu;
            minimizer = new Min_inv(f, gradf, mu);
        } else if(user_choice == 2){
            std::cout << "Please insert the value for sigma (suggested value between 0 and 0.5) : ";
            std::cin >> sigma;
            minimizer = new Min_arm(f, gradf, sigma);
        } else if(user_choice == 3) {
            std::cout << "Please insert the value for eta (suggested value = 0.9) : ";
            std::cin >> eta;
            std::cout << "Please insert the value for mu (suggested value = 0.2) : ";
            std::cin >> mu;
            minimizer = new Min_hvy(f, gradf, eta, mu);
        }
    }

    // call the solve method
    minimizer->solve();

    // print the requested output
    std::cout << "Number of iterations = " << minimizer->get_iter() << std::endl;

    std::cout << "Minimum value x = ";
    minimizer->get_x().print();

    std::cout << "f(x) = " << f(minimizer->get_x()) << std::endl;

    std::cout << "grad_f(x) = ";
    gradf(minimizer->get_x()).print();

    delete minimizer;

    return 0;
}
