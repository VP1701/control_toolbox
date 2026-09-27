// kalman_filter.cpp

#include "kalman_filter.hpp"


KalmanFilter::KalmanFilter(const matrix& A, const matrix& B, const matrix& C,
                           const matrix& Q, const matrix& R, const matrix& P, const matrix& x_hat) 
                           : nx(A.rows), ny(C.rows), nu(B.columns), A(A), B(B), C(C), Q(Q), R(R),
                            P(P), I(mops.eye(nx)), x_hat(x_hat) {
                                status_ = init();
                            }


    Status KalmanFilter::init() {
        if (A.rows != A.columns || A.rows != B.rows || C.columns != nx ||
            x_hat.rows != nx || P.rows != P.columns || P.rows != nx ||
            Q.rows != nx || R.rows != ny || x_hat.columns != 1) {
            return Status::DimensionMismatch;
        }

        return Status::Ok;
    }
   

    Status KalmanFilter::step(const matrix& z, const matrix& u, matrix& y_hat) {
        Status s = predict(u);
        if (s != Status::Ok) {
            return s;
        }
        return update(z, y_hat);
    }

    Status KalmanFilter::predict(const matrix& u) {
        if (status_ != Status::Ok) {
            return Status::NotInitialized;
        }
        if (u.rows != nu) {
            return Status::DimensionMismatch;
        }
        x_hat = A * x_hat + B * u;
        P = (A * P * A.T() + Q);
        return Status::Ok;
    }

    Status KalmanFilter::update(const matrix& z, matrix& y_hat) {
        if (status_ != Status::Ok) {
            return Status::NotInitialized;
        }
        if (z.rows != ny || y_hat.rows != ny || z.columns != 1 ||
            y_hat.columns != 1) {
            return Status::DimensionMismatch;
        }
        matrix innovation = z - C * x_hat;
        matrix S = C * P * C.T() + R;
        matrix K = P * C.T() * mops.inverse(S);
        x_hat = x_hat + K * innovation;
        P = (I - K * C) * P; // TODO: make more robust update
        y_hat = z - C * x_hat;
        return Status::Ok;
        
    }