// simplex.h

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <iostream>
#include <matrix.h>
#include <vector>

class simplex {
    private:
        int constraints;
        int variables;

        std::vector<int> basis;

        matrix A;
        matrix b;
        matrix c_trans;

        Matrix mops;

    public:
        Simplex(const matrix& A_in, const matrix& b_in, const matrix& c_trans_in)


}