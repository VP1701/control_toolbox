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
    matrix M(n, n);
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

matrix Matrix::lower_toeplitz_I(int T_size, int block_size) {
    

    matrix I = eye(block_size);
    matrix T(T_size,T_size);
    
    // make C_uId into toeplitz matrix which is a lower block triangular matrix where blocks are I
    for (int i = 0; i < T_size; i += block_size) {
        for (int j = (i / block_size) * block_size; j >= 0; j -= block_size) {
            T.set_block(i, j, I);
        }
    }

    return T;
}

matrix Matrix::create_block_diagonal(const matrix& A, int n_blocks) {
    int rows = A.rows * n_blocks;
    int columns = A.columns * n_blocks;
    matrix result(rows, columns);

    for (int i = 0; i < n_blocks; i++) {
        result.set_block(i * A.rows, i*A.columns, A);
    }
    return result;
}

matrix Matrix::diag(const matrix& A) {
    // returns diagonal elements of a square matrix
    int rows = A.rows;
    int columns = A.columns;
    matrix diag(rows,1); // columns vector to store diagonal elements

    for(int i = 0; i < A.rows; ++i) {
        diag(i,0) = A.data[i * columns + i];
    }
    return diag;
}


matrix Matrix::zeros(int r, int c) {
    matrix m(r, c);

    for (int i = 0; i < m.rows * m.columns; i++) {
        m.data[i] = 0.0;
    }

    return m;
};

matrix Matrix::ones(int r, int c) {
    matrix m(r, c);

    for (int i = 0; i < m.rows * m.columns; i++) {
        m.data[i] = 1.0;
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

// TODO: update these to be cache friendly later. doesnt matter for small matrices
matrix Matrix::multiply(const matrix& A, const matrix& B) {
    if (A.columns != B.rows) {
            throw std::runtime_error("Dimension don't match for multiplication");
    }
    matrix C(A.rows, B.columns);
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
    if (A.columns != B.columns || A.rows != B.rows) {
            throw std::runtime_error("Dimension don't match for addition");
    }
    matrix C(A.rows, A.columns);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.columns; j++) {
            C(i,j) = A(i,j) + B(i,j);
        }
    }
    return C;
}

matrix Matrix::subtraction(const matrix& A, const matrix& B) {
    if (A.columns != B.columns || A.rows != B.rows) {
            throw std::runtime_error("Dimension don't match for subtraction");
    }
    matrix C(A.rows, A.columns);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.columns; j++) {
            C(i,j) = A(i,j) - B(i,j);
        }
    }
    return C;
}

matrix Matrix::get_column(const matrix& A, int n) {
    int nx = A.rows;
    matrix column = zeros(nx, 1);
    for (int i = 0; i < nx; ++i) {
        column(i,0) = A(i,n);
    }
    return column;
}

matrix Matrix::get_row(const matrix& A, int n) {
    // wip
    int nx = A.rows;
    matrix column = zeros(nx, 1);
    
    return column;
}

matrix Matrix::lin_solve(const matrix& A, const matrix& b) {
    // solve x from A*x=b using cholesky decomposition
    matrix L = A.chol();
    int rows = L.rows;
    int columns = L.columns;
    matrix y(rows, 1);
    matrix x(rows, 1);

    // solve intermediary y
    for (int i = 0; i < rows; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L.data[i * columns + j]*y.data[j];
        }
        y.data[i] = (b.data[i] - sum) / L.data[i * columns + i];
    }

    // solve x without transposing L 
    // L_ij^T=L_ji
    for (int i = rows - 1; i >= 0; --i) { // start from final row
        double sum = 0.0;
        for (int j = i + 1; j < rows; ++j) {
            sum += L.data[j * columns + i]*x.data[j];
        }
        x.data[i] = (y.data[i] - sum) / L.data[i * columns + i];
    }
    return x;
}

matrix Matrix::lin_solve_LU(const matrix& A, const matrix& b) {
    int n = b.rows;
    matrix b_perm = b;

    matrix x(n,1);
    matrix y(n,1);
    
    LUResult result = A.LU();
    matrix& LU_mat = result.LU;
    int rows = LU_mat.rows;
    std::vector<int> permutations = result.permutations;

    // swap rows of b to match LU decomposition permutations
    for (int i = 0;i < n; ++i) {
        b_perm(i, 0) = b(permutations[i], 0);
    }
    double sum;
    // solve L*y = b
    for (int i = 0; i < n; ++i) {
        sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += y(j, 0) * LU_mat(i, j);
        }
        y(i, 0) = b_perm(i, 0) - sum;
        
    }

    // solve U*x = y
    for (int i = n - 1; i >= 0; --i) {
        sum = 0.0;
        for (int j = n - 1; j > i; --j) {
            sum += LU_mat(i, j) * x(j, 0);
        }
        x(i, 0) = ( y(i, 0) - sum ) / LU_mat(i, i);
    }

    return x;
}

matrix Matrix::lin_solve_chol(const matrix& L, const matrix& b) {
    // solve x from A*x=b using precomputed cholesky decomposition L
    
    int rows = L.rows;
    int columns = L.columns;
    matrix y(rows, 1);
    matrix x(rows, 1);

    // solve intermediary y
    for (int i = 0; i < rows; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L.data[i * columns + j]*y.data[j];
        }
        y.data[i] = (b.data[i] - sum) / L.data[i * columns + i];
    }

    // solve x without transposing L 
    // L_ij^T=L_ji
    for (int i = rows - 1; i >= 0; --i) { // start from final row
        double sum = 0.0;
        for (int j = i + 1; j < rows; ++j) {
            sum += L.data[j * columns + i]*x.data[j];
        }
        x.data[i] = (y.data[i] - sum) / L.data[i * columns + i];
    }
    return x;
}

matrix Matrix::chol_rank1_update(const matrix& L_M, const matrix& z, double beta) {
    // based on algortihm 3.1 from https://christian-igel.github.io/paper/AMERCMAUfES.pdf
    int n = L_M.rows;
    matrix result(n, L_M.columns);
    matrix omega = z;
    double b = 1.0;
    for (int i = 0; i < n; i++) {
        double l_ii = L_M(i,i);
        double l_ii_pow = l_ii*l_ii;
        double omega_pow = omega.data[i]*omega.data[i];
        double l_new_ii = std::sqrt(l_ii_pow+beta/b*omega_pow);
        double gamma = l_ii_pow*b+beta*omega_pow;
        for (int j = i + 1; j < n; j++) {
            double l_ji = L_M(j,i);
            omega.data[j] = omega.data[j] - omega.data[i]/l_ii*l_ji;
            double l_ji_new = (l_new_ii/l_ii)*l_ji+(l_new_ii*beta*omega.data[i])/gamma*omega.data[j];
            result(j,i) = l_ji_new;
        }
        b = b + beta * (omega.data[i]*omega.data[i])/l_ii_pow;
        result(i,i) = l_new_ii;
        
    }
    return result;
}

matrix Matrix::backslash(const matrix& A, const matrix& B) {
    // solves A * X = B with chol. only works for symmetric positive definite matrix A 
    // will update to use LU decomposition later to make function more general.

    matrix L = A.chol();
    matrix X(B.rows, B.columns);

    // split the problem to solve columns by columns so lin_solve_chol can be used

    for  (int i = 0; i < B.columns; ++i) {
        matrix b_col = get_column(B, i);
        
        matrix x_col = lin_solve_chol(L, b_col);

        X.set_block(0, i, x_col);
    }
    return X;
}

matrix Matrix::backslash_chol(const matrix& L, const matrix& B) {
    // solves A * X = B with chol. only works for symmetric positive definite matrix A.

    matrix X(B.rows, B.columns);

    // split the problem to solve columns by columns so lin_solve_chol can be used

    for  (int i = 0; i < B.columns; ++i) {
        matrix b_col = get_column(B, i);
        
        matrix x_col = lin_solve_chol(L, b_col);

        X.set_block(0, i, x_col);
    }
    return X;
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

        augmented.swap_rows(p, max_val_row);

        double pivot = augmented(p,p);
        augmented.multiply_row(p, 1.0 / pivot);

        // ELiminate column p from other rows
        for (int i = 0; i < n; ++i) {
            if (i != p) {
                double factor = augmented(i, p);

                augmented.add_multiple_of_row(i, p, -factor);
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

