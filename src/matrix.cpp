// matrix.cpp
#include "matrix.h"

matrix Matrix::swap_rows(matrix& A, int row1, int row2) {
    if (row1 == row2) return A;
    for (int j = 0; j < A.columns; ++j) {
        std::swap(A(row1,j), A(row2,j));
    }
    return A;
}

matrix Matrix::multiply_row(matrix& A, int row, double scalar) {
    if (scalar == 0.0) {
        std::cout << "Error: multiplying row by zero!" << "\n";
        return A;
    }
    for (int j = 0; j < A.columns; ++j) {
        A(row,j) *= scalar;
    }
    return A;
}

matrix Matrix::add_multiple_of_row(matrix& A, int dest, int src, double scalar) {
    for (int j = 0; j < A.columns; ++j) {
        A(dest, j) += scalar * A(src, j);
    }
    return A;
}

matrix Matrix::eye(int n) {
    matrix M;
    M.rows = n;
    M.columns = n;
    M.data = new double[M.rows * M.columns];
    for (int i = 0; i < M.rows; ++i) {
        for (int j = 0; j < M.columns; ++j) {
            if (i==j) {
                M(i,j)=1.0;
            } else{
                M(i,j)=0.0;
            }
        }
    }
    return M;
}

matrix Matrix::zeros(int r, int c) {
    matrix m;
    m.rows = r;
    m.columns = c;
    m.data = new double[m.rows * m.columns];

    for (int i = 0; i < m.rows * m.columns; i++) {
        m.data[i] = 0.0;
    }

    return m;
};

void Matrix::print(const matrix& A)  {
    for (int i = 0; i < A.rows; i++)  {
        for (int j = 0; j < A.columns; j++) {
            std::cout << A(i,j) << " ";
        }
        std::cout << "\n";
    }
}

matrix Matrix::multiply(const matrix& A, const matrix& B) {
    matrix C = zeros(A.rows, B.columns);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < B.columns; j++)  {
            for (int k = 0; k < A.columns; k++) {
                C(i,j) += A(i,k) * B(k, j);
            }
        }
    }
    return C;
}


matrix Matrix::addition(const matrix& A, const matrix& B) {
    matrix C = zeros(A.rows, A.columns);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.columns; j++) {
            C(i,j) = A(i,j) + B(i,j);
        }
    }
    return C;
}

matrix Matrix::subtraction(const matrix& A, const matrix& B) {
    matrix C = zeros(A.rows, A.columns);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.columns; j++) {
            C(i,j) = A(i,j) - B(i,j);
        }
    }
    return C;
}

matrix Matrix::get_column(const matrix& A, int n) {
    nx = A.rows;
    matrix column = zeros(nx, 1);
    for (int i = 0; i < nx; ++i) {
        column(i,0) = A(i,n);
    }
    return column;
}



matrix Matrix::inverse(const matrix& A) {
    int n = A.rows;
    if (A.rows != A.columns) {
        std::cout << "Not a square matrix, cannot invert." << "\n";
    }

    matrix augmented = zeros(n, 2*n);

    // copy A into the left part of the augmented matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented(i,j) = A(i,j);
        }
    }

    // make the right part of augmented matrix the identity matrix
    for (int i= 0; i < n; ++i) {
        augmented(i, i + n) = 1.0;
    }

    // Gauss-Jordan elimination

    for (int p = 0; p < n; ++p) {

        int max_val_row = p;
        for (int i = p + 1; i < n; ++i) {
            if (std::abs(augmented(i,p)) > std::abs(augmented(max_val_row, p))) {
                max_val_row = i;
            }
        }

        // Row swap to bring pivot to (p,p)
        if (std::abs(augmented(max_val_row, p)) < 1e-12) {
            std::cout << "Error: MAtrix is singular --> can not invert. \n";
            return A;
        }

        augmented = swap_rows(augmented, p, max_val_row);

        double pivot = augmented(p,p);
        augmented = multiply_row(augmented, p, 1.0 / pivot);

        // ELiminate column p from other rows
        for (int i = 0; i < n; ++i) {
            if (i != p) {
                double factor = augmented(i, p);

                add_multiple_of_row(augmented, i, p, -factor);
            }
        }
    }

    matrix inverted_A = zeros(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverted_A(i,j) = augmented(i,j+n);
        }
    }

    // implement with gauys-Jordan elimination
    return inverted_A;
}
