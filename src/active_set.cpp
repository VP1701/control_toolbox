// active_set.cpp

#include "active_set.h"


matrix Active_Set::get_active_rows(const matrix& A, const std::vector<int>& ws) {
    // returns a submatrix A_W that includes the rows of A that belong to the active set
    int rows = ws.size(); // amount of active constraint
    int columns = A.columns; 
    matrix A_W(rows, columns);

    for (int i = 0; i < rows;i++) {
        matrix A_r = A.get_row(ws[i]);
        A_W.set_block(i, 0, A_r);
    }

    return A_W;
};

void Active_Set::initialize_QP(const matrix& A, const matrix& b) {
    // the first feasible solution can be found using the simplex method (fletcher s 162)
    int m = A.rows;
    int n = A.columns;
    matrix A_simp(m, 2*m+n);
    matrix ct(2*n+m,1);
    
    matrix I = mops.eye(m);
    A_simp.set_block(0,0,A);
    A_simp.set_block(0,n,A*(-1.0));
    A_simp.set_block(0,2*n,I);

    // fill in ones to the end
    for (int i = 2*n; i < 2*n + m; i++) {
        ct.data[i] = 1.0;
    }

    constraint_types = std::vector<ConstraintType>(m, LEQ);

    solver = new Simplex(A_simp, b, ct, constraint_types);

    solver->solve();

    matrix x_simp = solver->return_solution();

    // check if x actually feasible. sum(r_i)=0 must be true
    // elements of r are the m last elements of x

    //extract r
    matrix r(1,m);
    for (int i = 0; i < m; i++) {
        r.data[i] = x_simp.data[2*n+i];
    }

    if (r.L1() > tol) {
        std::cout << "No feasible starting point for the QP. Reformulate!" << "\n";
        return;
    } 

    // extract x
    this->x = matrix(n,1);
    for (int i = 0; i < n; i++) {
        this->x.data[i] = x_simp.data[i] - x_simp.data[i+n];
    }

    // cosntruct initial active set 
    // a_t*x_1-b_i = 0 -> active
    // a_t*x_1-b_i > 0 -> inactive
    // a_t*x_1-b_i < 0 -> infeasible 

};

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