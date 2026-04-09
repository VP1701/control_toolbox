// active_set.cpp

#include "active_set.h"

void Active_Set::solve_kkt(const matrix& L_Q, const matrix& g, const matrix& A) {
    /*
    solves KKT system using the Schur complement method

    [Q A^T  [p      = [-g
     A O  ]  lamda]     0]

    */


    // Q * p + A^T * lamda = -g
    //  p + Q^-1 * A^T * lamda = -Q^-1 * g
    //  p = -Q^-1 * g - Q^-1 * A^T * lamda
    // solve v=Q^-1 * g from Q*v = g
    matrix v = mops.lin_solve_chol(L_Q, g);

    //  p = -v - Q^-1 * A^T * lamda

    // next we build M = A * Q^-1 * A^T which is the schur complement of G
    // we use a temproray matrix F = Q^-1 * A^T to avoin calculating matrix inverse
    matrix A_T = A.T();
    matrix F(A_T.rows,A_T.columns);
    F = mops.backslash_chol(L_Q, A_T);

    // M = A * F
    matrix M = A * F;

    // now we can write p = -v - F * lamda and substitute it into
    // A * p = 0 -> A * (-v - F * lamda) = 0
    // A * F * lamda = -A*v
    // M * lamda = -A*v

    // solve lamda from M * lamda = -A * v
    matrix L_M = M.chol();
    matrix mA_v = (A * v) * (-1.0); // -A * v
    lamda = mops.lin_solve_chol(L_M, mA_v);

    // solve p from Q * p = -(g + A^T * lamda)
    matrix rhs_p = (g + A.T() * lamda) * (-1.0);
    p = mops.lin_solve_chol(L_Q, rhs_p);

}

void Active_Set::solve_QP(const matrix& Q, const matrix& c, const matrix& A_ineq, const matrix& b_ineq) {
    
}
Active_Set::Active_Set(const matrix& Q, const matrix& c, const matrix& A_ineq, const matrix& b_ineq) {
    n = Q.rows;
    x = mops.zeros(n, 1); // solution vector

    L_Q = Q.chol();
}