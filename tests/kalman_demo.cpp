// kalman_demo.cpp

#include "matrix.h"
#include "kalman_filter.hpp"

#include <fstream>
#include <iostream>
#include <random>

/*
A = [1 1; 0 1];
B = [0.5; 1];
H = [1 0];
Q = 0.01*eye(2);
R = 0.1;
x_est = [0; 0];   % initial state
P = eye(2);
*/
Matrix mops;

int main() {
    const int nx = 2;
    const int nu = 1;
    const int ny = 1;

    matrix A(nx, nx);
    A(0, 0) = 1.0; A(0, 1) = 1.0;
    A(1, 0) = 0.0; A(1, 1) = 1.0;

    matrix B(nx, nu);
    B(0, 0) = 0.5;
    B(1, 0) = 1.0;

    matrix C(ny, nx);
    C(0, 0) = 1.0;
    C(0, 1) = 0.0;

    matrix Q = mops.eye(nx) * 0.01;

    matrix R(ny, ny);
    R(0, 0) = 0.1;

    matrix x0(nx, 1);
    x0(0, 0) = 0.0;
    x0(1, 0) = 0.0;

    matrix P0 = mops.eye(nx);

    KalmanFilter filter(A, B, C, Q, R, P0, x0);

    if (filter.status() != Status::Ok) {
        std::cerr << "filter init failed\n";
        return 1;
    }

    std::mt19937 rng(42);                          // fixed seed = reproducible
    std::normal_distribution<double> w(0.0, 0.1);  // process noise
    std::normal_distribution<double> v(0.0, 0.32); // measurement noise, sqrt(R)

    matrix x_true(nx, 1);   // starts at zero
    matrix u(nu, 1);
    u(0, 0) = 1.0;          // constant input

    matrix z(ny, 1);
    matrix y_hat(ny, 1);

    std::ofstream csv("kalman_results.csv");
    csv << "k,x_true_0,x_true_1,x_est_0,x_est_1,z,P00,P11\n";

    for (int k = 0; k < 100; ++k) {
        x_true = A * x_true + B * u;
        x_true(0, 0) += w(rng);
        x_true(1, 0) += w(rng);

        z = C * x_true;
        z(0, 0) += v(rng);

        if (filter.step(z, u, y_hat) != Status::Ok) {
            std::cerr << "filter step failed at k = " << k << "\n";
            return 1;
        }

        const matrix& x_est = filter.state();
        const matrix& P = filter.covariance();

        csv << k << ","
            << x_true(0, 0) << "," << x_true(1, 0) << ","
            << x_est(0, 0) << "," << x_est(1, 0) << ","
            << z(0, 0) << ","
            << P(0, 0) << "," << P(1, 1) << "\n";
    }

    return 0;
}