// kalman_filter.hpp

#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include "matrix.h"
#include "status.hpp"


class KalmanFilter {
    private:
        mutable Matrix mops;
        int nx;
        int ny;
        int nu;

        Status status_ = Status::NotInitialized; 

        matrix A;
        matrix B;
        matrix C;
        matrix Q;
        matrix R;
        matrix P; 
        matrix I;

        matrix x_hat;

        Status init();
    
    public:

        // Getters
        matrix state() const { return x_hat; }
        matrix covariance() const { return P; }
        Status status() const { return status_; }
        
        Status predict(const matrix& u);
        Status update(const matrix& z, matrix& y_hat);

        Status step(const matrix& z, const matrix& u, matrix& y_hat);
        KalmanFilter(const matrix& A, const matrix& B, const matrix& C,
                         const matrix& Q, const matrix& R, const matrix& P, const matrix& x_hat);

        ~KalmanFilter() = default;
        
};


#endif // kalman_filter.hpp