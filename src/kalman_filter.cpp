// kalman_filter.cpp

#include "kalman_filter.hpp"


KalmanFilter::KalmanFilter(const matrix& A, const matrix& B, const matrix& C,
                           const matrix& Q, const matrix& R) 
                           : nx(A.rows), A(A), B(B), C(C), Q(Q), R(R),
                            P(matrix(nx, nx)), I(mops.eye(nx)), x_hat(matrix(nx, 1)) {}


    matrix KalmanFilter::step(const matrix& z, const matrix& u) {
        x_hat = A * x_hat + B * u;
        P = (A * P * A.T() + Q);
        matrix innovation = z - C * x_hat;
        matrix S = C * P * C.T() + R;
        matrix K = P * C.T() * mops.inverse(S);
        x_hat = x_hat + K * innovation;
        P = (I - K * C) * P;
        matrix y_hat = z - C * x_hat;
        return y_hat;
    }