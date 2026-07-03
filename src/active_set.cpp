// active_set.cpp

#include "active_set.h"


matrix Active_Set::get_active_rows(const matrix& A, const std::set<int>& ws) {
    // returns a submatrix A_W that includes the rows of A that belong to the active set
    int rows = ws.size(); // amount of active constraint
    int columns = A.columns; 
    matrix A_W(rows, columns);


    int i = 0;
    for (int idx : ws) {
        matrix A_r = A.get_row(idx);
        A_W.set_block(i, 0, A_r);
        i++;
    }

    return A_W;
};

void Active_Set::initialize_QP() {
    // the first feasible solution can be found using the simplex method (fletcher s 162)
    delete solver;
    solver = nullptr;

    active_set.clear();

    int m = A.rows;
    int n = A.columns;
    matrix A_simp(m, 2*n+m);
    matrix ct(2*n+m,1);
    
    matrix I = mops.eye(m);
    A_simp.set_block(0,0,A*(-1.0));
    A_simp.set_block(0,n,A);
    A_simp.set_block(0,2*n,I);

    // fill in ones to the end
    for (int i = 2*n; i < 2*n + m; i++) {
        ct.data[i] = 1.0;
    }

    constraint_types = std::vector<ConstraintType>(m, LEQ);

    solver = new Simplex(A_simp, b*(-1.0), ct, constraint_types);

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
        //std::cout << "No feasible starting point for the QP. Reformulate!" << "\n";
        return;
    } 

    // extract x
    this->x = matrix(n,1);
    for (int i = 0; i < n; i++) {
        this->x.data[i] = x_simp.data[i] - x_simp.data[i+n];
    }

    // cosntruct initial active set 
    // a_t*x-b_i = 0 -> active
    // a_t*x-b_i > 0 -> inactive
    // a_t*x-b_i < 0 -> infeasible 

    for (int i = 0; i < m; i++) {
        matrix A_r = A.get_row(i);
        w_residual = (A_r * x).scalar() - b.data[i];

        if (std::abs(w_residual) < tol) {
            active_set.insert(i);
        } else if (w_residual < -tol) {
            //std::cout << "The problem is infeasible" << "\n";
            return;
        }
    }

    if (!active_set.empty()) {
        matrix A_W = get_active_rows(A, active_set);
        matrix F = mops.backslash_chol(L_Q, A_W.T());
        matrix M = A_W * F;
        L_M = M.chol();
    } else {
        L_M = matrix(0, 0);
    }

    std::cout << "QP initialized \n";
    
};

matrix Active_Set::return_solution() {
    return x;
}

void Active_Set::solve_kkt(const matrix& g, const matrix& A) {
    /*
    solves KKT system using the Schur complement method

    [Q A^T  [p      = [-g
     A O  ]  lamda]     0]

    */

    //std::cout << "Solving KKT" << "\n";
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
    // colesky rank 1 update makes this more efficient
    //matrix M = M + a_i * Q_inv * a_i.T();
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


// helper for debug. move to be part of matrix struct
bool has_nan(const matrix& A) {
    for (int i = 0; i < A.rows*A.columns; ++i)
        if (std::isnan(A.data[i])) return true;
    return false;
}

void Active_Set::solve() {
    
    for (int i = 0; i < max_iters; i++) {
        matrix A_W = get_active_rows(A, active_set);

        bool is_nan = has_nan(A_W);
        // (is_nan) std::cout << "A_W has nan\n";
        matrix g = Q * x + c;
        if (active_set.empty()){
            p = mops.lin_solve_chol(L_Q, g*(-1.0));
        } else {
            solve_kkt(g, A_W);
        }
        
        

        // update x

        // check if ||p|| = 0

        if (p.L2() < tol) {
            if (active_set.empty()) return;
            // all lamda must be non-negative
            bool not_optimal = false;
            for (int k = 0; k < lamda.rows; k++) {
                if (lamda.data[k] >= tol) {
                    not_optimal = true;
                }
            }

            if (not_optimal) {
                // blands rule for active set. pick the lowest index -> never cycle
                double largest = -1e100;
                int index = -1;
                for (int a = 0; a < lamda.rows * lamda.columns; a++) {
                    double val = lamda.data[a];
                    if (val > tol) {
                        index = a;
                        break;
                    }
                }

                auto it = active_set.begin();
                std::advance(it, index);
                
                active_set.erase(*it);
                /*if (!active_set.empty()) {
                    matrix A_W = get_active_rows(A, active_set);
                    matrix F = mops.backslash_chol(L_Q, A_W.T());
                    L_M = (A_W * F).chol();
                } else {
                    L_M = matrix(0,0);
                }*/
                continue;
            } else {
                return;
            }
        }

        // calculate step size alpha
        double result = 1e9;
        int index = -1;
        for (int j = 0; j < n_const; j++) {
            if (active_set.count(j) == 0) {
                matrix A_r = A.get_row(j);
                double den = (A_r*p).scalar();
                if (den < -tol) {
                    double num = b.data[j] - (A_r*x).scalar();
                
                    double frac = num/den;
                    if (frac < result) {
                        result = frac;
                        index = j;
                    }
                }
            }
        }

        double alpha = std::min(1.0, result);

        if (alpha < 1.0 && index != -1) {
            /*
            // rank1 update to L_M when adding constraints
            matrix a_j = A.get_row(index).T();
            z = Q_inv * a_j;
            //Golub & Van Loan Section 6.5.4
            L_M  = mops.chol_rank1_update(L_M, z, 1.0);*/
            active_set.insert(index);
        }

        x = x + p * alpha;
    }

    //std::cout << "Maximum number of iteraions reached. Solution not optimal!" << "\n";
}

void Active_Set::reset(const matrix& b, const matrix& c_new) {
   
    this->b = b;
    this->c = c_new;
    active_set.clear();
    for (int i = 0; i < A.rows; i++) {
        double res = (A.get_row(i) * x - b.get_row(i)).scalar();
        if ( res < -tol) {
            initialize_QP();
            return;
        }
        if (std::abs(res) < tol) {
            active_set.insert(i);
        }
    }

    if (!active_set.empty()) {
        matrix A_W = get_active_rows(A, active_set);
        matrix F = mops.backslash_chol(L_Q, A_W.T());
        matrix M = A_W * F;
        L_M = M.chol();
    } else {
        L_M = matrix(0, 0);
    }

}

Active_Set::Active_Set(const matrix& Q, const matrix& c, const matrix& A_ineq, const matrix& b_ineq) {
    solver = nullptr;
    n = Q.rows;
    x = mops.zeros(n, 1); // solution vector
    this->Q = Q;
    Q_inv = mops.inverse(Q);
    A = A_ineq;
    b = b_ineq;
    n_const = A.rows;
    this->c = c;
    L_Q = Q.chol();

    initialize_QP();

}

Active_Set::~Active_Set() {
    delete solver;
}