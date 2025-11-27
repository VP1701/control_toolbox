// simplex.cpp


#include <iostream>
#include <simplex.h>


matrix Simplex::get_basis() const {

}


Simplex::simplex(const matrix& A_in, const matrix& b_in, const matrix& c_trans_in) {
    

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
            A(i,i) = A_in(i,j);
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
    c_trans_big = zeros(n_big, 1);

    // copy A into left side of A_big
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A_big(i,i) = A(i,j);
        }
    }

}