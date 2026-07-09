// matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath> // for sqrt
#include <vector>
#include <stdexcept>

struct LUResult;

struct matrix {
    int rows;
    int columns;
    double* data;

    matrix() : rows(0), columns(0), data(nullptr) {}

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

    matrix& operator=(matrix&& other) noexcept {
        if (this == &other) return *this;

        delete[] data;

        rows = other.rows;
        columns = other.columns;
        data = other.data;

        other.data = nullptr;

        other.rows = 0;
        other.columns = 0;

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

    double L1() const {
        double sum = 0.0;
        int iters = rows * columns;

        for (int i = 0; i < iters; i++) {
            sum += std::abs(data[i]);
        }
        return sum;
    }

    double L2() const {
        double sum = 0.0;
        int iters = rows * columns;

        for (int i = 0; i < iters; i++) {
            sum += std::pow(data[i],2);
        }
        return std::sqrt(sum);
    }

    double scalar() const {
        if (rows != 1 || columns != 1) {
            std::cout << "Not a 1x1 matrix. Cant convert to scalar!" << "\n";
        }

        return data[0];
    }

    double min() const {
        // Returns smalles t element from matrix
        double smallest = 1e100;
        for (int i = 0; i < rows * columns; i++) {
            double val = data[i];
            if (val < smallest) {
                smallest = val;
            }
        }
        return smallest;
    }

    matrix chol() const {
        matrix L(rows, columns);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L.data[i * columns + k] * L.data[j * columns + k]; //L[i][k] * L[j][k]
                }

                if (i == j) {
                    L(i,j) = std::sqrt(data[i * columns + i] - sum);//L[i][j] = sqrt(A[i][i] - sum)
                } else {
                    L(i,j) = (data[i * columns + j] - sum) / (L.data[j * columns + j]); //L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum))
                }
            }
        }
        return L;
    }

    matrix get_row(int n) const {
        // return the n:th row of the matrix
        matrix A_r(1, columns);
        // check feasibility
        if (n >= rows || n < 0) {
            std::cout << "ERROR: row index out of bounds, can't extract!" << "\n";
            return A_r;
        }
        //extract row
        for (int i = 0; i < columns; i++) {
            A_r.data[i] = data[n * columns + i];
        }

        return A_r;
    }

    matrix& swap_rows(int row_ind1, int row_ind2) {
        // Check that indices are different
        if (row_ind1 == row_ind2) return *this;
        
        for (int i = 0; i < columns; ++i) {
            std::swap(data[row_ind1 * columns + i], data[row_ind2 * columns + i]);
        }
        return *this;
    }

    matrix& multiply_row(int row_ind, double scalar) {
        if (scalar == 0.0) {
        std::cout << "Error: multiplying row by zero!" << "\n";
        }
        for (int i = 0; i < columns; ++i) {
            data[row_ind * columns + i] = data[row_ind * columns + i] * scalar;
        }
        return *this;
    }

    matrix& add_multiple_of_row(int dest_row, int src_row, double scalar) {
        if (scalar == 0.0) return *this;
        
        for (int i = 0; i < columns; ++i) {
            data[dest_row * columns + i] += data[src_row * columns + i] * scalar;
        }

        return *this;
    }

    // set block of matrix to some matrix
    matrix& set_block(int start_row, int start_column, const matrix& A) {
        if (start_row < 0 || start_column < 0 || start_row + A.rows > rows || start_column + A.columns > columns) {
            std::cout << "Impossible substitution! Check your dimesnsion!" << "\n";
            return *this;
        }
        for (int i = 0; i < A.rows; ++i) {
            for (int j = 0; j < A.columns; ++j) {

                data[start_row * columns + start_column + i * columns + j] = A.data[i * A.columns + j];
            }
        }
        return *this;
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
    
    // trick to make scalar multiplication work. but need to always m,ultiply from the right :(

    matrix operator * (double scalar) const {
        matrix result(rows, columns);
        for (int i = 0; i < rows * columns; ++i) {
            result.data[i] = data[i] * scalar;
        }
        return result;
    }

    matrix operator + (const matrix& B) {
        matrix result(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result.data[i * columns + j] = data[i * columns + j] + B.data[i * columns + j];
            }
        }
        return result;
    }

    matrix operator - (const matrix& B) {
        matrix result(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result.data[i * columns + j] = data[i * columns + j] - B.data[i * columns + j];
            }
        }
        return result;
    }

};

struct LUResult {
    matrix LU;
    std::vector<int> permutations;
};

class Matrix {
    public:
        // elementary row operations
        matrix swap_rows(matrix& A, int row1, int row2);
        matrix multiply_row(matrix& A, int row, double scalar);
        matrix add_multiple_of_row(matrix& A, int dest, int src, double scalar);

        // functions for creating matrices
        matrix lower_toeplitz_I(int T_size, int block_size);
        matrix create_block_diagonal(const matrix& A, int n_blocks);
        matrix diag(const matrix& A);
        matrix eye(int n);
        matrix zeros(int r, int c);
        matrix ones(int r, int c);

        matrix backslash(const matrix& A, const matrix& B);
        matrix backslash_chol(const matrix& L, const matrix& B);
        matrix lin_solve(const matrix& A, const matrix& b);
        matrix lin_solve_LU(const matrix& A, const matrix& b);
        matrix lin_solve_chol(const matrix& L, const matrix& b);
        matrix chol_rank1_update(const matrix& L_M, const matrix& z, double beta);
        
        void print(const matrix& A);
        matrix multiply(const matrix& A, const matrix& B);
        matrix addition(const matrix& A, const matrix& B);
        matrix subtraction(const matrix& A, const matrix& B);
        matrix inverse(const matrix& A);
        
        matrix get_column(const matrix& A, int n);
        matrix get_row(const matrix& A, int n);
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

inline LUResult matrix::LU() const {
    matrix LU = *this;

    if (rows != columns) {
        throw std::runtime_error("Matrix must be square to LU decompose.");
        }

    std::vector<int> permutations(rows); 
    // Fill permutations vector with row indices
    for (int k = 0; k < rows; ++k) {
        permutations[k] = k;
    }
    // Find the largest (absolute) value in a column.
    
    for (int i = 0; i < columns; ++i) {
        int max_val_row = i;
        double max_val = std::abs(LU(max_val_row, i));
        for ( int j = i + 1; j < rows; ++j) {
            double val = std::abs(LU(j, i));
            if (val > max_val) {
                max_val = val;
                max_val_row = j;
            }
        }
        // Check for singularity. prvent division by zero
        if (max_val < 1e-12) {
            throw std::runtime_error("Matrix is singular, cannot LU decompose.");
        }
        // Swap largest value to the diagonal by doing a row permutation
        LU.swap_rows(max_val_row, i);
        std::swap(permutations[max_val_row], permutations[i]);

        // zero out everything below diagonal
        // go through all of the rows below diagonal 
        // and multiply them by m times diagonal row
        for (int j = i + 1;j < rows; ++j) {
            // LU(j,i) - m * L(i,i) = 0 --> m = LU(j,i) / L(i,i)
            double m = LU(j, i) / LU(i, i);
            LU.add_multiple_of_row(j, i, -m);
            LU(j, i) = m;
        }
    }

    return {LU, permutations};
}

#endif // MATRIX_H