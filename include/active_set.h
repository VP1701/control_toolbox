// active_set.h

#ifndef ACTIVE_SET_H
#define ACTIVE_SET_H

#include <iostream>
#include <chrono>
#include <matrix.h>
#include "simplex.h"
#include <vector>
#include <set>

class Active_Set {
    private:
        mutable Matrix mops;
        matrix Q;
        matrix Q_inv;
        matrix A;
        matrix c;
        matrix b;

        int n_const;

        // KKT
        matrix p; // search direction
        matrix lamda; // lagrangian multipliers

        int n; // amount of variables
        matrix x;
        matrix L_Q;
        matrix L_M;
        matrix z;

        // residual for checking if constraint belongs to active set
        double w_residual = 0.0;

        std::set<int> active_set; // active set indices

        matrix get_active_rows(const matrix& A, const std::set<int>& ws);


        //helpers for simplex algorithm used inside
        std::vector<ConstraintType> constraint_types;
        Simplex* solver;

        double tol = 1e-6;
        double convergence_tol = 1e-6;
        int max_iters = 500;

    public:
        void solve_kkt(const matrix& g, const matrix& A);
        void initialize_QP();
        void solve();
        void reset(const matrix& b, const matrix& c_new);
        matrix return_solution();
        Active_Set(const matrix& Q, const matrix& c, const matrix& A_ineq, const matrix& b_ineq);
        ~Active_Set();


    };

#endif // ACTIVE_SET_H