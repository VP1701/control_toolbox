// kalman_filter.hpp

#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include "matrix.h"

class KalmanFilter {
    private:
        mutable Matrix mops;
        int nx;
        matrix A;
        matrix B;
        matrix C;
        matrix Q;
        matrix R;
        matrix P; 
        matrix I;

        matrix x_hat;
    
    public:

        matrix step(const matrix& z, const matrix& u);
        KalmanFilter(const matrix& A, const matrix& B, const matrix& C,
                         const matrix& Q, const matrix& R);

        ~KalmanFilter() = default;
        
};


#endif // kalman_filter.hpp