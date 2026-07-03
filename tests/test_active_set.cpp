// test_active_set.cpp

#include <iostream>
#include "matrix.h"
#include "active_set.h"

void run_test(const char* name, const matrix& Q, const matrix& c, const matrix& A, const matrix& b) {
    std::cout << "\n";
    std::cout << "==================================================\n";
    std::cout << "TEST: " << name << "\n";
    std::cout << "==================================================\n";

    Active_Set solver(Q, c, A, b);
    solver.solve();

    matrix x = solver.return_solution();

    std::cout << "Solution:\n";
    for (int i = 0; i < x.rows; ++i) {
        std::cout << "x_" << i << " = " << x(i, 0) << "\n";
    }

    matrix Qx = Q * x;
    double quad_term = (x.T() * Qx).scalar();
    double lin_term  = (c.T() * x).scalar();
    double obj = 0.5 * quad_term + lin_term;

    std::cout << "Objective value (0.5*x'Qx + c'x): " << obj << "\n";
}

int main() {
    Matrix mops;

    // ------------------------------------------------------------
    // Test 1: unconstrained optimum already feasible (no active constraints)
    // min  x1^2 + x2^2 - 2x1 - 2x2   ->  Q = 2I, c = [-2, -2]
    // Active_Set expects constraints in the form A x >= b, so
    // "x1 + x2 <= 10" becomes "-x1 - x2 >= -10", and x1>=0, x2>=0 stay as-is.
    // Unconstrained minimum is at x = (1, 1), well inside the feasible region.
    // ------------------------------------------------------------
    {
        matrix Q = mops.zeros(2, 2);
        Q(0,0) = 2; Q(1,1) = 2;

        matrix c = mops.zeros(2, 1);
        c(0,0) = -2; c(1,0) = -2;

        matrix A = mops.zeros(3, 2);
        A(0,0) = -1; A(0,1) = -1;
        A(1,0) = 1;  A(1,1) = 0;
        A(2,0) = 0;  A(2,1) = 1;

        matrix b = mops.zeros(3, 1);
        b(0,0) = -10; b(1,0) = 0; b(2,0) = 0;

        run_test("Interior optimum, no active constraints", Q, c, A, b);
        std::cout << "expected: x = [1, 1], opt_cost = -2\n";
    }

    // ------------------------------------------------------------
    // Test 2: Q = I, single active constraint (simple projection)
    // min 0.5*(x1^2+x2^2) - x1 - x2  -> Q = I, c = [-1,-1]
    // "x1 + x2 <= 1" becomes "-x1 - x2 >= -1"; x1>=0, x2>=0 stay as-is.
    // Unconstrained minimum (1,1) violates x1+x2<=1, so it becomes active.
    // Projecting (1,1) onto the line x1+x2=1 gives (0.5, 0.5).
    // ------------------------------------------------------------
    {
        matrix Q = mops.zeros(2, 2);
        Q(0,0) = 1; Q(1,1) = 1;

        matrix c = mops.zeros(2, 1);
        c(0,0) = -1; c(1,0) = -1;

        matrix A = mops.zeros(3, 2);
        A(0,0) = -1; A(0,1) = -1;
        A(1,0) = 1;  A(1,1) = 0;
        A(2,0) = 0;  A(2,1) = 1;

        matrix b = mops.zeros(3, 1);
        b(0,0) = -1; b(1,0) = 0; b(2,0) = 0;

        run_test("Single active constraint (projection)", Q, c, A, b);
        std::cout << "expected: x = [0.5, 0.5], opt_cost = -0.75\n";
    }

    // ------------------------------------------------------------
    // Test 3: classic Nocedal & Wright example (Numerical Optimization)
    // min (x1-1)^2 + (x2-2.5)^2   ->  Q = 2I, c = [-2, -5]  (+const 7.25, dropped)
    // Original (textbook) form is <=; rewritten here in the A x >= b form
    // that Active_Set expects, by negating each row:
    //  x1-2x2  >= -2   (was -x1+2x2 <= 2)
    // -x1-2x2  >= -6   (was  x1+2x2 <= 6)
    // -x1+2x2  >= -2   (was  x1-2x2 <= 2)
    //  x1      >=  0
    //      x2  >=  0
    // ------------------------------------------------------------
    {
        matrix Q = mops.zeros(2, 2);
        Q(0,0) = 2; Q(1,1) = 2;

        matrix c = mops.zeros(2, 1);
        c(0,0) = -2; c(1,0) = -5;

        matrix A = mops.zeros(5, 2);
        A(0,0) = 1;  A(0,1) = -2;
        A(1,0) = -1; A(1,1) = -2;
        A(2,0) = -1; A(2,1) = 2;
        A(3,0) = 1;  A(3,1) = 0;
        A(4,0) = 0;  A(4,1) = 1;

        matrix b = mops.zeros(5, 1);
        b(0,0) = -2; b(1,0) = -6; b(2,0) = -2; b(3,0) = 0; b(4,0) = 0;

        run_test("Nocedal & Wright classic QP", Q, c, A, b);
        std::cout << "expected: x = [1.4, 1.7], f(x) = (x1-1)^2+(x2-2.5)^2 = 0.8\n";
        std::cout << "(internal objective without the +7.25 constant is -6.45)\n";
    }

    // ------------------------------------------------------------
    // Test 4: 3-variable QP with a non-diagonal Q block
    // min 0.5*x'Qx + c'x,  Q = [[4,1,0],[1,2,0],[0,0,2]], c = [-6,-3,-4]
    // "x1+x2+x3 <= 3" becomes "-x1-x2-x3 >= -3"; x1,x2,x3 >= 0 stay as-is.
    // Unconstrained optimum sums to ~4.14, so the sum constraint is active
    // (only that one -- nonnegativity stays slack at the solution).
    // Expected values below confirmed via compare_active_set.m (MATLAB quadprog).
    // ------------------------------------------------------------
    {
        matrix Q = mops.zeros(3, 3);
        Q(0,0) = 4; Q(0,1) = 1; Q(0,2) = 0;
        Q(1,0) = 1; Q(1,1) = 2; Q(1,2) = 0;
        Q(2,0) = 0; Q(2,1) = 0; Q(2,2) = 2;

        matrix c = mops.zeros(3, 1);
        c(0,0) = -6; c(1,0) = -3; c(2,0) = -4;

        matrix A = mops.zeros(4, 3);
        A(0,0) = -1; A(0,1) = -1; A(0,2) = -1;
        A(1,0) = 1;  A(1,1) = 0;  A(1,2) = 0;
        A(2,0) = 0;  A(2,1) = 1;  A(2,2) = 0;
        A(3,0) = 0;  A(3,1) = 0;  A(3,2) = 1;

        matrix b = mops.zeros(4, 1);
        b(0,0) = -3; b(1,0) = 0; b(2,0) = 0; b(3,0) = 0;

        run_test("3D QP with coupled variables", Q, c, A, b);
        std::cout << "expected (MATLAB quadprog): x = [1.1333, 0.4, 1.4667], opt_cost = -8.5333\n";
    }

    std::cout << "\nALL TESTS COMPLETED!\n";
    return 0;
}