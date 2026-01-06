// matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

struct matrix {
    int rows;
    int columns;
    double* data;

    // constructor
    matrix(int r, int c) : rows(r), columns(c) {
        data = new double[r * c]();
    }

    // destructor
    ~matrix() {
        delete[] data;
    }

    // copy contructor
    matrix(const matrix& other) : rows(other.rows), columns(other.columns) {
        data = new double[rows * columns];
        for (int i = 0; i < rows * columns; ++i) {
            data[i] = other.data[i];
        }
    }

    // copy matrix operator
    matrix& operator=(const matrix& other) {
        // check if copying self
        if (this == &other) return *this;

        delete[] data;

        // copy
        rows = other.rows;
        columns = other.columns;
        data = new double[rows * columns];
        for (int i = 0; i < rows* columns; ++i) {
            data[i] = other.data[i];
        }

        return *this;
    }

    // matrix transpose
    matrix transpose() const {
        matrix result(columns, rows);

        for (int i = 0; i < rows; ++i) {
            int offset = i * columns;
            for (int j = 0; j < columns; ++j) {
                result.data[j * rows + i] = data[offset + j];
            }
    }

    return result;
    }
    double& operator()(int r, int c) {
        return data[r * columns + c];
    }

    double operator()(int r, int c) const {
        return data[r * columns + c];
    }
};

class Matrix {
    public:
        // elementary row operations
        matrix swap_rows(matrix& A, int row1, int row2);
        matrix multiply_row(matrix& A, int row, double scalar);
        matrix add_multiple_of_row(matrix& A, int dest, int src, double scalar);




        matrix eye(int n);
        matrix zeros(int r, int c);
        void print(const matrix& A);
        matrix multiply(const matrix& A, const matrix& B);
        matrix addition(const matrix& A, const matrix& B);
        matrix subtraction(const matrix& A, const matrix& B);
        matrix inverse(const matrix& A);
        matrix get_column(const matrix& A, int n);
};

inline matrix operator*(const matrix& A, const matrix& B) {
    Matrix mops;
    return mops.multiply(A, B);
}

inline matrix operator+(const matrix& A, const matrix& B) {
    Matrix mops;
    return mops.addition(A, B);
}

inline matrix operator-(const matrix& A, const matrix& B) {
    Matrix mops;
    return mops.subtraction(A, B);
}
#endif // MATRIX_H