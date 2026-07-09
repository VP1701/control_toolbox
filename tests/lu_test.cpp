#include <iostream>
#include <cmath>
#include "matrix.h"

// Prints A*x and compares it against b, reporting max absolute residual.
void check_solution(const char* name, const matrix& A, const matrix& x, const matrix& b) {
    Matrix mops;
    matrix Ax = mops.multiply(A, x);
    double max_resid = 0.0;
    for (int i = 0; i < b.rows; ++i) {
        max_resid = std::max(max_resid, std::abs(Ax(i, 0) - b(i, 0)));
    }
    std::cout << "TEST: " << name << "\n";
    std::cout << "  max |A*x - b| = " << max_resid << "  -> "
               << (max_resid < 1e-9 ? "PASS" : "FAIL") << "\n";
}

int main() {
    Matrix mops;

    {
        // Classic partial-pivoting example: requires two row swaps to factor.
        matrix A = mops.zeros(3, 3);
        A(0,0)=2; A(0,1)=1; A(0,2)=1;
        A(1,0)=4; A(1,1)=3; A(1,2)=3;
        A(2,0)=8; A(2,1)=7; A(2,2)=9;

        matrix b = mops.zeros(3, 1);
        b(0,0)=5; b(1,0)=10; b(2,0)=20;

        matrix x = mops.lin_solve_LU(A, b);
        check_solution("3x3 requiring pivoting", A, x, b);
    }

    {
        // Identity matrix: LU should need no elimination at all, x should equal b exactly.
        matrix A = mops.eye(4);
        matrix b = mops.zeros(4, 1);
        b(0,0)=1; b(1,0)=2; b(2,0)=3; b(3,0)=4;

        matrix x = mops.lin_solve_LU(A, b);
        check_solution("identity matrix", A, x, b);
    }

    {
        // Diagonal matrix: no pivoting needed, but exercises non-trivial divisions.
        matrix A = mops.zeros(3, 3);
        A(0,0)=2; A(1,1)=5; A(2,2)=10;
        matrix b = mops.zeros(3, 1);
        b(0,0)=4; b(1,0)=15; b(2,0)=40;

        matrix x = mops.lin_solve_LU(A, b);
        check_solution("diagonal matrix", A, x, b);
    }

    {
        // Larger system with no special structure, cross-checked against Matrix::inverse
        // (a completely independent solve path) rather than a hand-picked b.
        matrix A = mops.zeros(5, 5);
        double vals[5][5] = {
            { 6,  1, -2,  0,  3},
            { 2,  8,  1, -1,  0},
            {-1,  0,  7,  2,  1},
            { 3, -2,  0,  9, -1},
            { 0,  1,  2, -3,  5}
        };
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                A(i,j) = vals[i][j];

        matrix b = mops.zeros(5, 1);
        double bvals[5] = {1, 2, 3, 4, 5};
        for (int i = 0; i < 5; ++i) b(i,0) = bvals[i];

        matrix x_lu = mops.lin_solve_LU(A, b);
        matrix x_inv = mops.multiply(mops.inverse(A), b);

        check_solution("5x5 vs Matrix::inverse (lin_solve_LU)", A, x_lu, b);
        check_solution("5x5 vs Matrix::inverse (inverse path)", A, x_inv, b);

        double max_diff = 0.0;
        for (int i = 0; i < 5; ++i)
            max_diff = std::max(max_diff, std::abs(x_lu(i,0) - x_inv(i,0)));
        std::cout << "  max |x_lu - x_inv| = " << max_diff << "  -> "
                   << (max_diff < 1e-6 ? "PASS" : "FAIL") << "\n";
    }

    {
        // Singular matrix: LU() must throw, not silently return garbage.
        matrix A = mops.zeros(2, 2);
        A(0,0)=1; A(0,1)=2;
        A(1,0)=2; A(1,1)=4;   // row 2 = 2 * row 1 -> singular

        std::cout << "TEST: singular matrix throws\n";
        try {
            LUResult result = A.LU();
            std::cout << "  FAIL: no exception thrown\n";
        } catch (const std::runtime_error& e) {
            std::cout << "  PASS: caught exception: " << e.what() << "\n";
        }
    }

    {
        // Non-square matrix: LU() must throw, not read/write out of bounds.
        matrix A = mops.zeros(2, 3);

        std::cout << "TEST: non-square matrix throws\n";
        try {
            LUResult result = A.LU();
            std::cout << "  FAIL: no exception thrown\n";
        } catch (const std::runtime_error& e) {
            std::cout << "  PASS: caught exception: " << e.what() << "\n";
        }
    }

    return 0;
}