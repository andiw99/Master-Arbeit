//
// Created by andi on 20.07.23.
//

#ifndef LEARNINGPROJECT_QUENCHCORRLENGTH_H
#define LEARNINGPROJECT_QUENCHCORRLENGTH_H

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

void printMatrixXd(const Eigen::MatrixXd& matrix) {
    std::cout << "Matrix values:\n";
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            std::cout << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

#endif //LEARNINGPROJECT_QUENCHCORRLENGTH_H
