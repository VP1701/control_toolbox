// matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

struct matrix {
    int rows;
    int columns;
    double* data;

    double& operator()(int r, int c) {
        return data[r * columns + c];
    }

    double operator()(int r, int c) const {
        return data[r * columns + c];
    }
};

class Matrix {
    public:

        matrix zeros(int r, int c);
        void print(matrix A);
        matrix multiply(matrix A, matrix B);
        matrix addition(matrix A, matrix B);
};

inline matrix operator*(const matrix& A, const matrix& B) {
    Matrix mops;
    return mops.multiply(A, B);
}

inline matrix operator+(const matrix& A, const matrix& B) {
    Matrix mops;
    return mops.addition(A, B);
}
#endif // MATRIX_H