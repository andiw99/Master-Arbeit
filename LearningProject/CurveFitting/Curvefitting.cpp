//
// Created by andi on 14.07.23.
//

#include "Curvefitting.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include <iostream>
#include "../Header/Helpfunctions and Classes.h"

using namespace std;

template<typename scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
    typedef scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

};

struct ExponentialDecayFunctor : Functor<double> {
private:
    // Observations for a sample.
    Eigen::MatrixXd C_x;
public:
    ExponentialDecayFunctor(const Eigen::MatrixXd &C_x): Functor(C_x.cols(), C_x.rows()), C_x(C_x) {
        cout << "rows: " << C_x.rows() << endl;
        cout << "cols: " << C_x.cols() << endl;
    }

    int operator()(const Eigen::VectorXd &paras, Eigen::VectorXd &fvec) const {
        // fvec is the vector with residuals
        for(int i = 0; i < this->values(); i++) {
            fvec(i) = C_x(i, 1) - paras(0) * exp(- C_x(i, 0) / paras(1));
        }
        return 0;
    }
};

struct LorentzianPeakFunctor : Functor<double> {
private:
    // Observations for a sample.
    Eigen::MatrixXd X_Y;
public:
    LorentzianPeakFunctor(const Eigen::MatrixXd &X_Y): Functor(X_Y.cols(), X_Y.rows()), X_Y(X_Y) {
        cout << "rows: " << X_Y.rows() << endl;
        cout << "cols: " << X_Y.cols() << endl;
    }

    int operator()(const Eigen::VectorXd &paras, Eigen::VectorXd &fvec) const {
        // fvec is the vector with residuals
        for(int i = 0; i < this->values(); i++) {
            fvec(i) = X_Y(i, 1) - paras(0) * 2 * paras(1) /
                    (1 + 4 * M_PI * M_PI * X_Y(i, 0) * X_Y(i,0) * paras(1) * paras(1));
        }
        return 0;
    }
};

struct DivergenceRightFunctor : Functor<double> {
private:
    // Observations for a sample.
    Eigen::MatrixXd X_Y;
public:
    DivergenceRightFunctor(const Eigen::MatrixXd &X_Y): Functor(X_Y.cols(), X_Y.rows()), X_Y(X_Y) {
        cout << "rows: " << X_Y.rows() << endl;
        cout << "cols: " << X_Y.cols() << endl;
    }

    int operator()(const Eigen::VectorXd &paras, Eigen::VectorXd &fvec) const {
        // fvec is the vector with residuals

        // paras(0) = T-C
        // paras(1) = xi0
        // paras(2) = nu

        // writing the functor that optimizes to the right of the divergence, meaning we have negative
        // reduced temperatures
        for(int i = 0; i < this->values(); i++) {
            // double eps = (paras(0) - X_Y(i, 0)) / paras(0);
            // will this work like this?
            fvec(i) = X_Y(i, 1) - paras(1) / pow(-(paras(0) - X_Y(i, 0)) / paras(0) , paras(2));
        }
        return 0;
    }
};

struct NumericalExpDecayDiffFunctor : Eigen::NumericalDiff<ExponentialDecayFunctor> {
public:
    NumericalExpDecayDiffFunctor(const Eigen::MatrixXd &C_x) : Eigen::NumericalDiff<ExponentialDecayFunctor>(C_x)
    {}
};

class Fitter {
private:
    int kNumObservations;
public:
    Fitter() {

    }

    template <class Residual>
    vector<double> fit() {

    }
};

template <template<class valuetype> class vec, class value_type>
vec<value_type> average_two(const vec<value_type> &a, const vec<value_type> &b) {
    // vec needs to be initializable like a vector
    vec<value_type> result(a.size());

    // averaging
    for(int i = 0; i < a.size(); i++) {
        result[i] = (a[i] + b[i]) / (value_type) 2.0;
    }

    return result;
}

template <class vec>
vec average_two(const vec &a, const vec &b) {
    // vec needs to be initializable like a vector
    vec result(a.size());

    // averaging
    for(int i = 0; i < a.size(); i++) {
        result[i] = (a[i] + b[i]) / 2.0;
    }

    return result;
}

template <template<typename valuetype> typename vec, typename value_type>
Eigen::MatrixXd construct_matrix(const vec<value_type> &a, const vec<value_type> &b) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd>(a.data(), a.size());
    Eigen::VectorXd b_eigen = Eigen::Map<Eigen::VectorXd>(b.data(), b.size());

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

template <typename vec>
Eigen::MatrixXd construct_matrix(const vec &a, const vec &b) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a.data(), a.size());
    Eigen::VectorXd b_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), b.size());

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

Eigen::MatrixXd construct_matrix(vector<double> &a, vector<double> &b) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a.data(), (long)a.size());
    Eigen::VectorXd b_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), (long)b.size());

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

Eigen::MatrixXd construct_matrix(Eigen::VectorXd &a, Eigen::VectorXd &b) {

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a, b;

    return result_matrix;
}

void printMatrixXd(const Eigen::MatrixXd& matrix) {
    std::cout << "Matrix values:\n";
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            std::cout << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

#include <random>

// Generate value pairs for a Lorentzian-type peak with noise
std::vector<double> generateLorentzianData(std::vector<double> &xValues)
{
    std::vector<double> yValues(xValues.size());

    // Parameters for the Lorentzian peak
    double xi = 5.0;  // Full width at half maximum (FWHM)
    double amplitude = 20.0;  // Peak amplitude

    // Random number generator for noise
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(0.0, 0.01);  // Gaussian noise distribution

    // Generate x-values and corresponding y-values with Lorentzian shape and added noise
    for (int i = 0; i < xValues.size(); ++i)
    {
        double x = xValues[i];  // You can modify the x-values based on your requirements

        // Lorentzian shape
        double y = amplitude * (2 * xi / (1 + 4 * M_PI * M_PI * x * x * xi * xi));

        // Add Gaussian noise
        double noise = distribution(generator);
        y += noise;

        // Store the value pairs
        xValues[i] = x;
        yValues[i] = y;
    }

    // Return the generated value pairs
    return yValues;
}

std::vector<double> generateDivergenceData(std::vector<double> &xValues)
{
    std::vector<double> yValues(xValues.size());

    // Parameters for the Lorentzian peak
    double amplitude = 5.0;  // Full width at half maximum (FWHM)
    double nu = 1.0;  // Peak amplitude
    double Tc = 100.0;

    // Random number generator for noise
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(0.0, 0.1);  // Gaussian noise distribution

    // Generate x-values and corresponding y-values with Lorentzian shape and added noise
    for (int i = 0; i < xValues.size(); ++i)
    {
        double x = xValues[i];  // You can modify the x-values based on your requirements

        double eps = (Tc - x) / Tc;
        cout << eps << endl;
        // Lorentzian shape
        double y = amplitude / (pow(abs(eps), nu));

        // Add Gaussian noise
        double noise = distribution(generator);
        // y += noise;

        // Store the value pairs
        xValues[i] = x;
        yValues[i] = y;
    }

    // Return the generated value pairs
    return yValues;
}

template<typename T>
std::vector<T> linspace(T start_in, T end_in, int num_in)
{

    std::vector<T> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back((T)(start + delta * i));
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}


int main() {
    cout << "hello world" << endl;
    // What is even the situation that we have?
    // we have had values for distances and the correlation function, but i don't think that we want to fit that now
    // but we can try it just for the sake of practice
    std::vector<double> xValues = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<double> yValues = {8.24, 6.92, 4.78, 3.51, 2.17, 1.67, 1.12, 0.85, 0.41, 0.19};

    Eigen::MatrixXd X_Y_vals(10, 2);
    X_Y_vals = construct_matrix(xValues, yValues);
    printMatrixXd(X_Y_vals);
    Eigen::VectorXd params(2);
    params << 1.0, 1.0;
    ExponentialDecayFunctor functor(X_Y_vals);
    cout << functor.values() << endl;
    cout << functor.inputs() << endl;
    Eigen::NumericalDiff<ExponentialDecayFunctor> numericalDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<ExponentialDecayFunctor>> lm(numericalDiff);

    Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(params);
    std::cout << "status: " << status << std::endl;
    std::cout << "amplitude: " << params(0) << std::endl;
    std::cout << "xi: " << params(1) << std::endl;


    xValues = linspace( -3.0, 3.0, 20);
    yValues = generateLorentzianData(xValues);
    X_Y_vals = construct_matrix(xValues, yValues);
    printMatrixXd(X_Y_vals);
    params << 1.0, 1.0;
    LorentzianPeakFunctor functor2(X_Y_vals);
    cout << functor2.values() << endl;
    cout << functor2.inputs() << endl;
    Eigen::NumericalDiff<LorentzianPeakFunctor> numericalDiff2(functor2);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LorentzianPeakFunctor>> lm2(numericalDiff2);

    Eigen::LevenbergMarquardtSpace::Status status2 = lm2.minimize(params);
    std::cout << "status: " << status2 << std::endl;
    std::cout << "amplitude: " << params(0) << std::endl;
    std::cout << "xi: " << params(1) << std::endl;


    xValues = linspace( 101.0,  200.0, 20);
    yValues = generateDivergenceData(xValues);
    X_Y_vals = construct_matrix(xValues, yValues);
    printMatrixXd(X_Y_vals);
    Eigen::VectorXd params2(3);
    params2 << 110.0, 5.0, 1.2;
    DivergenceRightFunctor functor3(X_Y_vals);
    cout << functor3.values() << endl;
    cout << functor3.inputs() << endl;
    Eigen::NumericalDiff<DivergenceRightFunctor> numericalDiff3(functor3);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<DivergenceRightFunctor>> lm3(numericalDiff3);
    lm3.parameters.maxfev = 10000000;
    Eigen::LevenbergMarquardtSpace::Status status3 = lm3.minimize(params2);
    std::cout << "status: " << status3 << std::endl;
    std::cout << "T_c: " << params2(0) << std::endl;
    std::cout << "xi0: " << params2(1) << std::endl;
    std::cout << "nu: " << params2(2) << std::endl;

    return 0;
}