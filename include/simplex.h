// simplex.h

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <iostream>
#include <matrix.h>
#include <vector>

class Simplex {
    private:
        int m; // amount of constraints
        int n; // amount of original svariables
        int n_big; // amount of variables in big-Ms

        std::vector<int> basis;


        matrix A;
        matrix b;
        matrix c_trans;
        matrix A_big;
        matrix c_trans_big;
        matrix x;
        std::vector<int> artificial_indices;

        double M = 1e9; // Big-M value


        Matrix mops;

        // functions

        matrix get_basis() const;
        calculate_current_solution() const;
        void print_solution() const;s

    public:
        // constructor
        Simplex(const matrix& A_in,
                const matrix& b_in,
                const matrix& c_trans_in);
        
        // Solves the optimization problem
        void solve();



}