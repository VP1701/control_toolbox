// active_set.h

#ifndef ACTIVE_SET_H
#define ACTIVE_SET_H

#include <iostream>
#include <chrono>
#include <matrix.h>
#include "simplex.h"
#include <vector>

class Active_Set {
    private:
        mutable Matrix mops;
        matrix A;
        matrix A_W;
        matrix p;
        matrix lamda;
        int n; // amount of variables
        matrix x;
        matrix L_Q;
        int iters = 100;

        std::vector<int> working_set; // active set indices

        matrix get_active_rows(const matrix& A, const std::vector<int>& ws);


        //helpers for simplex algorithm
        std::vector<ConstraintType> constraint_types;
        Simplex* solver;

        double tol = 1e-9;
        int max_iters = 100;

        matrix x_init;
        double mu = 1e6; // for L1 exact penalty coef

    public:
        void solve_kkt(const matrix& Q, const matrix& g, const matrix& A);
        void initialize_QP(const matrix& A, const matrix& b);
        void solve_QP(const matrix& Q, const matrix& c, const matrix& A_ineq, const matrix& b_ineq);
        Active_Set(const matrix& Q, const matrix& c, const matrix& A_ineq, const matrix& b_ineq);
        ~Active_Set();


    };

#endif // ACTIVE_SET_H