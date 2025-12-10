// main.cpp
// FULL VALIDATION SUITE FOR YOUR SIMPLEX SOLVER

#include <iostream>
#include <vector> // Required for std::vector
#include "matrix.h"
#include "simplex.h" // Assumes ConstraintType enum is here

// Define the run_test function to accept the vector of constraint types
void run_test(const char* name, const matrix& A, const matrix& b, const matrix& c, 
              const std::vector<ConstraintType>& types) { // ADDED ARGUMENT
    std::cout << "\n";
    std::cout << "==================================================\n";
    std::cout << "TEST: " << name << "\n";
    std::cout << "==================================================\n";

    // Call the updated constructor
    Simplex solver(A, b, c, types); 
    // solver.solve() is called in constructor

    
}

int main() {
    Matrix mops;

    {
        matrix A = mops.zeros(2,2);
        A(0,0)=2; A(0,1)=3;
        A(1,0)=-1; A(1,1)=1;

        matrix b = mops.zeros(2,1);
        b(0,0)=6; b(1,0)=1;

        matrix c = mops.zeros(1,2);
        c(0,0)=-1; c(0,1)=-3;

        // Constraint Types: Both are LEQ
        std::vector<ConstraintType> types = {LEQ, LEQ};

        run_test("lecture notes", A, b, c, types);
        std::cout << "answer should be: " << "x_transpose = [0.6, 1.6], opt_cost = -5.4" << "\n";
    }

    {
        matrix A_inf = mops.zeros(2,2);
        // Original: x₁ + x₂ ≥ 2 -> Used as -x₁ -x₂ ≤ -2 in the setup
        A_inf(0,0)=1; A_inf(0,1)=1; 
        A_inf(1,0)= 1; A_inf(1,1)= -1; 

        matrix b_inf = mops.zeros(2,1);
        b_inf(0,0)= 1;
        b_inf(1,0)= -2;

        matrix c_inf = mops.zeros(1,2);
        c_inf(0,0)=-1; c_inf(0,1)=-1; 
        std::vector<ConstraintType> types_inf = {LEQ, LEQ}; 
        run_test("INFEASIBLE problem", A_inf, b_inf, c_inf, types_inf); 
    }

    {
        matrix A = mops.zeros(1,2);
        A(0,0)= -1; A(0,1)= -1; 

        matrix b = mops.zeros(1,1);
        b(0,0)= -1;

        matrix c = mops.zeros(1,2);
        c(0,0)=-1; c(0,1)=-2;

        std::vector<ConstraintType> types = {LEQ}; 

        run_test("UNBOUNDED problem", A, b, c, types);
    }

    {
        matrix A = mops.zeros(2,2);
        A(0,0)=1; A(0,1)=4; 
        A(1,0)=1; A(1,1)=2; 

        matrix b = mops.zeros(2,1);
        b(0,0)=8; b(1,0)=4;

        matrix c = mops.zeros(1,2);
        c(0,0)=-3; c(0,1)=-9;

        std::vector<ConstraintType> types = {LEQ, LEQ};
        run_test("Degeneracy test", A, b, c, types);
    }

    {
    matrix A = mops.zeros(8,10);
    // Row-wise filling (values from the seeded rand above)
    double A_data[8][10] = {
        {1.484, 0.937, 2.744, 1.095, 0.632, 2.173, 1.918, 0.049, 2.501, 1.723},
        {0.294, 2.181, 1.639, 0.917, 2.884, 0.425, 1.273, 2.660, 0.167, 1.984},
        {2.463, 0.751, 1.852, 2.094, 0.388, 1.629, 2.987, 0.883, 1.526, 0.279},
        {1.172, 2.695, 0.604, 1.483, 2.350, 0.976, 0.251, 1.840, 2.712, 0.668},
        {0.895, 1.329, 2.108, 0.447, 1.763, 2.584, 0.732, 1.195, 0.339, 2.876},
        {2.017, 0.583, 1.964, 2.331, 0.108, 1.472, 2.695, 0.624, 1.883, 0.950},
        {1.706, 2.249, 0.417, 1.628, 2.973, 0.862, 1.354, 2.108, 0.775, 1.491},
        {0.628, 1.883, 2.540, 0.296, 1.017, 2.762, 0.194, 1.650, 2.429, 0.583}
    };

    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 10; ++j)
            A(i,j) = A_data[i][j];

    matrix b = mops.zeros(8,1);
    double b_vals[8] = {18.743, 19.294, 19.826, 18.472, 19.105, 19.683, 20.351, 19.628};
    for (int i = 0; i < 8; ++i) b(i,0) = b_vals[i];

    matrix c = mops.zeros(1,10);
    double c_vals[10] = {-2.351, -3.827, -1.194, -3.506, -2.683, -1.927, -3.148, -2.074, -1.583, -3.992};
    for (int j = 0; j < 10; ++j) c(0,j) = c_vals[j];

    std::vector<ConstraintType> types(8, LEQ);

    run_test("Large-ish feasible problem", A, b, c, types);
    }


    {
    matrix A = mops.zeros(6, 2);  // 6 rows, 2 variables
    A(0,0)=1;    A(0,1)=1;
    A(1,0)=1;    A(1,1)=0.25;
    A(2,0)=1;    A(2,1)=-1;
    A(3,0)=-0.25;A(3,1)=-1;
    A(4,0)=-1;   A(4,1)=-1;
    A(5,0)=-1;   A(5,1)=1;

    matrix b = mops.zeros(6,1);
    b(0,0)=2;
    b(1,0)=1;
    b(2,0)=2;
    b(3,0)=1;
    b(4,0)=-1;   // ← this is ≤ -1 → will be flipped to ≥ 1
    b(5,0)=2;

    matrix c = mops.zeros(1,2);
    c(0,0)=-1;
    c(0,1)=-1.0/3.0;

    // ALL ARE INEQUALITIES (≤)
    std::vector<ConstraintType> types = {LEQ, LEQ, LEQ, LEQ, LEQ, LEQ};

    //run_test("Klee-Minty cube (2D version)", A, b, c, types);
    }
    

    {
    Matrix mops;

    // 9 variables
    matrix A = mops.zeros(12, 9);
    matrix b = mops.zeros(12, 1);
    matrix c = mops.zeros(1, 9);

    std::vector<ConstraintType> types = {
        LEQ, LEQ, LEQ,
        LEQ, LEQ, LEQ,
        LEQ, LEQ, LEQ,
        LEQ, LEQ, LEQ
    };

    // ============================================================
    // OBJECTIVE: MAX profit into MIN (-profit)
    // ============================================================
    c(0,0) = -(125 - 58);  // x_LR
    c(0,1) = -(132 - 58);  // x_LM
    c(0,2) = -(140 - 58);  // x_LP

    c(0,3) = -(125 - 65);  // x_MR
    c(0,4) = -(132 - 65);  // x_MM
    c(0,5) = -(140 - 65);  // x_MP

    c(0,6) = -(125 - 74);  // x_HR
    c(0,7) = -(132 - 74);  // x_HM
    c(0,8) = -(140 - 74);  // x_HP

    // ============================================================
    // FEEDSTOCK constraint
    // ============================================================
    A(0,0)=1; A(0,1)=1; A(0,2)=1;   b(0,0)=20000;
    A(1,3)=1; A(1,4)=1; A(1,5)=1;   b(1,0)=20000;
    A(2,6)=1; A(2,7)=1; A(2,8)=1;   b(2,0)=10000;

    // ============================================================
    // MARKET DEMAND constraint
    // ============================================================
    // 4: -(x_LR + x_MR + x_HR) <= -20000
    A(3,0)=-1; A(3,3)=-1; A(3,6)=-1;  b(3,0)=-20000;

    // 5: -(x_LM + x_MM + x_HM) <= -10000
    A(4,1)=-1; A(4,4)=-1; A(4,7)=-1;  b(4,0)=-10000;

    // 6: -(x_LP + x_MP + x_HP) <= -5000
    A(5,2)=-1; A(5,5)=-1; A(5,8)=-1;  b(5,0)=-5000;

    // ============================================================
    // OCTANE constraint
    // ============================================================
    // 7: (3 x_LR - 2 x_MR - 6 x_HR) <= 0
    A(6,0)=3; A(6,3)=-2; A(6,6)=-6;  b(6,0)=0;

    // 8: (5 x_LM - 4 x_HM) <= 0
    A(7,1)=5; A(7,7)=-4;             b(7,0)=0;

    // 9: (9 x_LP + 4 x_MP) <= 0
    A(8,2)=9; A(8,4)=4;              b(8,0)=0;

    // ============================================================
    // SULFUR constraint
    // ============================================================
    A(9,0)=0.20; A(9,3)=-0.60; A(9,6)=-0.90;   b(9,0)=0;
    A(10,1)=0.50; A(10,4)=-0.30; A(10,7)=-0.60; b(10,0)=0;
    A(11,2)=0.80; A(11,8)=-0.30;                 b(11,0)=0;

    run_test("Refinery blending optimization (min version)", A, b, c, types);
    }
    std::cout << "\nALL TESTS COMPLETED!\n";

    return 0;
}