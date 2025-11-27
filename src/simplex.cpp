// simplex.cpp


#include <iostream>
#include <simplex.h>


matrix Simplex::get_basis() const {

    matrix B = mops.zeros(m, m);

    for (int i = 0; i < m; ++i) {
        int A_col = basis[i];
        for (j = 0; j < m; ++j) {
            B(j,i) = A(j, A_col);
        }
    }
    return B;
}

Simplex::calculate_current_solution() const {
    // calculate xB = B^(-1) * b

    matrix B = get_basis();
    matrix B_inv = mops.inverse(B);
    matrix xB = B_inv * b;

    matrix x = mops.zeros(n_big, 1);

    for (int i = 0; i < m; ++i) {
        int idx = basis[i];
        x(idx, 0) = xB(i, 0);
    }

}

void Simplex::print_solution() const {

    std::cout << "Printing solution to the simplex" << "\n";
    mops.print(x);
    std::cout << "Optimal value: " << "\n";
    matrix opt = c_trans * x;
    mops.print(opt);
}

void Simplex::solve() {
    // todo
    std::cout << "SOLVER NOT IMPLEMENTED YET" << "\n";
}

Simplex::Simplex(const matrix& A_in, const matrix& b_in, const matrix& c_trans_in) {
    

    int m = A_in.rows;
    int n = A_in.columns;

    // initialize matrices for optimization
    matrix A = mops.zeros(m,n);
    matrix b = mops.zeros(m,1);
    matrix c_trans = mops.zeros(1,n);

    // check that inputs have valid dimensions
    if (b_in.rows != m) {
        std::cout << "Invalid size of A or b. different amount of rows" << "\n";
    }

    if (b_in.columns != 1) {
        std::cout << "Invalid amount of columns on b. must have one column!" << "\n";
    }

    if (c_trans_in.rows != 1) {
        std::cout << "Invalid amount of rows in c_trans! Must have one row!" << "\n";
    }

    if (c_trans_in.columns != n) {
        std::cout << "Invalid amount of columns in A or c_trans! Must ahve same amount of columns!" << "\n";
    }

    // copy data from input matrices
    // REDO after implementing matrix struct in a smarter way

    // copy A_in data to A
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i,j) = A_in(i,j);
        }
    }

    // copy b_in data to b
    for (int i = 0; i < m; ++i) {
        b(i,0) = b_in(i,0);
    }

    // copy c_trans_in data to c_trans
    for (int i = 0; i < m; ++i) {
        c_trans(0,i) = c_trans_in(0,i);
    }

    // construct big-M matrices
    int n_big = n + m; 
    A_big = zeros(m, n_big);
    c_trans_big = zeros(1, n_big);

    // copy A into left side of A_big
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A_big(i,j) = A(i,j);
        }
    }

    // identity to right side of A_big
    artificial_indices.clear();
    for  (int i = 0; i < m; ++i) {
        A_big(i, n + i) = 1.0;
        int col = n + i;
        c_trans_big(0, col) = M;
        artificial_indices.push_back(col);
    }

    for (int j = 0; j < n; ++j) {
        c_trans_big(0, j) = c_trans_in(0, j);
    }

    basis.resize(m);
    for (int i = 0; i < m; ++i) {
        basis[i] = n + i;
    }

    A = A_big;
    c_trans = c_trans_big;

    std::cout << "Simplex initialized" << "\n";s

    calculate_current_solution();
    print_solution();
    
}