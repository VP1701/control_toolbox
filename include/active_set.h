// active_set.h

#ifndef ACTIVE_SET_H
#define ACTIVE_SET_H

#include <iostream>
#include <chrono>
#include <matrix.h>
#include <vector>

class Active_Set {
    private:
        mutable Matrix mops;
        matrix A;
        matrix p;
        matrix lamda;
        int n; // amount of variables
        matrix x;
        matrix L_Q;
        int iters = 100;

        double tol = 1e-9;
        int max_iters = 100;

        matrix x_init;
        double mu = 1e6; // for L1 exact penalty coef

    public:
        void solve_kkt(const matrix& Q, const matrix& g, const matrix& A);
        void solve_QP(const matrix& Q, const matrix& c, const matrix& A_ineq, const matrix& b_ineq);
        Active_Set(const matrix& Q, const matrix& c, const matrix& A_ineq, const matrix& b_ineq);
};

#endif // ACTIVE_SET_H