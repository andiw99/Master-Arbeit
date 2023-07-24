//
// Created by weitze73 on 08.06.23.
//

#include "EigenCurveFitting.h"
#include <iostream>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>

// Define the functor for curve fitting
struct CurveFittingFunctor {
    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
        // Compute the residual vector
        fvec(0) = x(0) * std::exp(-x(1) * 0.1) - 5.0;
        fvec(1) = x(0) * std::exp(-x(1) * 0.2) - 3.0;
        fvec(2) = x(0) * std::exp(-x(1) * 0.3) - 1.0;
        return 0;
    }

    int values() const {
        return 3;  // Number of residuals
    }


    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const {
        // Compute the Jacobian matrix
        fjac(0, 0) = std::exp(-x(1) * 0.1);
        fjac(0, 1) = -0.1 * x(0) * std::exp(-x(1) * 0.1);
        fjac(1, 0) = std::exp(-x(1) * 0.2);
        fjac(1, 1) = -0.2 * x(0) * std::exp(-x(1) * 0.2);
        fjac(2, 0) = std::exp(-x(1) * 0.3);
        fjac(2, 1) = -0.3 * x(0) * std::exp(-x(1) * 0.3);
        return 0;
    }

};

int main() {
    // Create an instance of the functor
    CurveFittingFunctor functor;

    // Create initial parameter values
    Eigen::VectorXd x(2);
    x << 1.0, 0.1;

    // Create an instance of Levenberg-Marquardt optimization algorithm
    Eigen::LevenbergMarquardt<CurveFittingFunctor> lm(functor);

    // Set the parameters for the optimization algorithm
    lm.parameters.maxfev = 100;  // maximum number of function evaluations
    lm.parameters.xtol = 1e-6;   // tolerance for the solution

    // Run the optimization
    lm.minimize(x);

    // Print the optimized parameter values
    std::cout << "Optimized parameters: " << x.transpose() << std::endl;

    return 0;
}