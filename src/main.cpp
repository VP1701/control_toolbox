// main.cpp


#include <iostream>
#include <vector> 
#include "matrix.h"
#include "simplex.h"
#include "mpc.h"
#include "fstream"

Matrix mops;


int main() {
    

    // 1. Setup the Dense SPD Matrix A
    matrix A(3, 3);
    A(0,0) = 4.0;   A(0,1) = 12.0;  A(0,2) = -16.0;
    A(1,0) = 12.0;  A(1,1) = 37.0;  A(1,2) = -43.0;
    A(2,0) = -16.0; A(2,1) = -43.0; A(2,2) = 98.0;

    std::cout << "--- 0. Factorization Check ---" << "\n";
    matrix L = A.chol();
    std::cout << "Computed L (Should be lower triangular: 2, 6, 1, -8, 5, 3): \n";
    mops.print(L);

    // ==========================================
    // TEST 1 & 2: Vector Solvers (lin_solve)
    // ==========================================
    matrix b(3, 1);
    b(0,0) = -20.0;
    b(1,0) = -43.0;
    b(2,0) = 192.0;

    std::cout << "\n--- 1. Testing lin_solve(A, b) ---" << "\n";
    matrix x1 = mops.lin_solve(A, b);
    std::cout << "Expected: [1, 2, 3]^T \nGot: \n";
    mops.print(x1);

    std::cout << "\n--- 2. Testing lin_solve_chol(L, b) ---" << "\n";
    matrix x2 = mops.lin_solve_chol(L, b);
    std::cout << "Expected: [1, 2, 3]^T \nGot: \n";
    mops.print(x2);


    // ==========================================
    // TEST 3 & 4: Matrix Solvers (backslash)
    // ==========================================
    matrix B(3, 2);
    B(0,0) = -20.0; B(0,1) = 28.0;
    B(1,0) = -43.0; B(1,1) = 80.0;
    B(2,0) = 192.0; B(2,1) = -141.0;

    std::cout << "\n--- 3. Testing backslash(A, B) ---" << "\n";
    matrix X1 = mops.backslash(A, B);
    std::cout << "Expected: \n[1, 0]\n[2, 1]\n[3, -1]\nGot: \n";
    mops.print(X1);

    std::cout << "\n--- 4. Testing backslash_chol(L, B) ---" << "\n";
    matrix X2 = mops.backslash_chol(L, B);
    std::cout << "Expected: \n[1, 0]\n[2, 1]\n[3, -1]\nGot: \n";
    mops.print(X2);

    return 0;
}