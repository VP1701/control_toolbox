// simplex.h

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <iostream>
#include <matrix.h>
#include <vector>

class Simplex {
    private:
        int m; // amount of constraints
        int n; // amount of variables

        std::vector<int> basis;

        matrix A;
        matrix b;
        matrix c_trans;
        matrix A_big;
        matrix c_trans_big;
        std::vector<int> artificial_indices;

        double M = 1e9; // Big-M value

        Matrix mops;

        // functions

        matrix Simplex::get_basis() const;



    public:
        Simplex::simplex(const matrix& A_in, const matrix& b_in, const matrix& c_trans_in);

        void solve();


}