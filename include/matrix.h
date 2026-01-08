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
    matrix T() const {
        matrix result(columns, rows);

        for (int i = 0; i < rows; ++i) {
            int offset = i * columns;
            for (int j = 0; j < columns; ++j) {
                result.data[j * rows + i] = data[offset + j];
            }
        }
        return result;
    }

    // set block of matrix to some matrix
    void set_block(int start_row, int start_column, const matrix& A) {
        if (start_row < 0 || start_column < 0 || start_row + A.rows > rows || start_column + A.columns > columns) {
            std::cout << "Impossible substitution! Check your dimesnsion!" << "\n";
            return;
        }
        for (int i = 0; i < A.rows; ++i) {
            for (int j = 0; j < A.columns; ++j) {

                data[start_row * columns + start_column + i * columns + j] = A.data[i * A.columns + j];
            }
        }
    }

    double& operator()(int r, int c) {
        return data[r * columns + c];
    }

    double operator()(int r, int c) const {
        return data[r * columns + c];
    }

    matrix operator * (const matrix& B) {
        matrix result(rows, B.columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < B.columns; j++)  {
                for (int k = 0; k < columns; k++) {
                    result.data[i * result.columns + j] += data[i * columns + k] * B.data[k * B.columns + j];
                }
            }
        }
        return result;
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