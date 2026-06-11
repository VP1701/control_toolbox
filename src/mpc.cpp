// mpc.cpp
#include "mpc.h"

extern Matrix mops;

matrix MPC::construct_A_a(const matrix& A, const matrix& B) {
    // C_a = [A, B;O I]

    int columns = A.columns + B.columns;
    int rows = A.columns + B.columns;
    int c = A.columns;
    int r = rows - A.rows;
    matrix O = mops.zeros(r, c);
    matrix I = mops.eye(B.columns);
    matrix A_a(rows, columns);
    A_a.set_block(0, 0, A);
    A_a.set_block(0, c, B);
    A_a.set_block(A.rows, 0, O);
    A_a.set_block(A.rows, c, I);
    return A_a;
}

matrix MPC::construct_B_a(const matrix& A, const matrix& B) {
    // B_a = [B;I]
    int columns = B.columns;
    int rows = A.columns + B.columns;
    matrix I = mops.eye(columns);
    matrix B_a(rows, columns);
    B_a.set_block(0, 0, B);
    B_a.set_block(B.rows, 0, I);
    return B_a;
}

matrix MPC::construct_C_a(const matrix& C, const matrix& B) {
    // C_a = [C, O]
    int rows = C.rows;
    int columns = C.columns + B.columns;
    matrix O = mops.zeros(rows, B.columns);
    matrix C_a(rows, columns);
    C_a.set_block(0, 0, C);
    C_a.set_block(0, C.columns, O);
    return C_a;
}

matrix MPC::construct_P(const matrix& A,  const matrix& C, int horizon) {
    /*
    Creates matrix S_c for prediction the state values x(k+1,...,k+h_p) = S_x * x(k) + dS_u * u(k-1) + S_u * u(k,...,k+h_u-1)
    This is for the case where all states are measurable so C matrix is included in S_x which is not valid
    for cases where not all states can be measured.
    */
    int n = C.rows;
    int columns = A.columns;
    matrix S_x(horizon * n, columns);
    matrix current_A = C * A;
    S_x.set_block(0, 0, current_A);
    
    for (int i = 1;i < horizon; ++i) {
        current_A = current_A*A;
        S_x.set_block(i * n, 0, current_A);
    }

    return S_x;
}

matrix MPC::construct_P_x(const matrix& A, int horizon) {
    /*
    Creates matrix S_c for prediction the state values x(k+1,...,k+h_p) = S_x * x(k) + dS_u * u(k-1) + S_u * u(k,...,k+h_u-1)
    This is for the case where all states are measurable so C matrix is included in S_x which is not valid
    for cases where not all states cannot be measured.
    */
    int n = A.rows;
    int columns = A.columns;
    matrix P_x(horizon * n, columns);
    matrix current_A = A;
    P_x.set_block(0, 0, current_A);
    
    for (int i = 1;i < horizon; ++i) {
        current_A = current_A*A;
        P_x.set_block(i * n, 0, current_A);
    }

    return P_x;
}

matrix MPC::construct_H(const matrix& A,const matrix& B, const matrix& C, int horizon) {
    // 
    int n = C.rows;
    int m = B.columns;
    int rows = horizon * n;
    int columns = horizon * m;

    matrix H(rows, columns);
    

    for (int i = 0; i < rows; i += n) {
        //matrix cumsum = C * B;
        int A_power = 0;
        for (int j = (i / n) * m; j >= 0; j -= m) {
            A_power = (i / n) - (j / m);
            
            //mops.print(S_u);
            //std::cout << "j: " << j << "\n";
            matrix CA = C; 
            for (int k = 0; k < A_power; k++) {
                //std::cout << "k: " << k << "\n";
                CA = CA * A;
                //std::cout << "CA: " << "\n";
                //mops.print(CA);
            }
            
            
            matrix cumsum = CA * B;
            //std::cout << "cumsum: " << "\n";
            //mops.print(cumsum);

            H.set_block(i, j, cumsum);
        }
    }
    return H;
}

matrix MPC::construct_H_x(const matrix& A,const matrix& B, int horizon) {
    int n = A.rows;
    int m = B.columns;
    int rows = horizon * n;
    int columns = horizon * m;

    matrix H(rows, columns);

    for (int i = 0; i < rows; i += n) {
        matrix cumsum = B;
        for (int j = (i / n) * m; j >= 0; j -= m) {
            //std::cout << "j: " << j << "\n";
            H.set_block(i, j, cumsum);
            //mops.print(S_u);
            cumsum = A * cumsum;
        }
    }
    return H;
}

matrix MPC::construct_constraint_C_du() {
    matrix I = mops.eye(nu*h);
    matrix C_du(2*h*nu,nu*h);
    matrix minusI = I*(-1);
    C_du.set_block(0,0,I);
    C_du.set_block(h*nu,0,minusI);
    return C_du;
}

matrix MPC::construct_constraint_C_u() {
    int rows = h*nu;

    matrix I = mops.eye(nu);
    matrix C_uId(rows,rows);
    matrix C_u(2*rows,rows);
    
    // make C_uId into toeplitz matrix which is a lower block triangular matrix where blocks are I
    for (int i = 0; i < rows; i += nu) {
        for (int j = (i / nu) * nu; j >= 0; j -= nu) {
            C_uId.set_block(i, j, I);
        }
    }
    C_u.set_block(0,0,C_uId);
    C_u.set_block(rows,0,C_uId*(-1));
    return C_u;
    
}

matrix MPC::construct_constraint_C_y() {
    int columns = H.columns;
    int rows = H.rows;
    matrix C_y(2 * rows, columns);

    C_y.set_block(0,0,H);
    C_y.set_block(rows,0,H*(-1));
    return C_y;
}

matrix MPC::construct_constraint_A(const matrix& C_du, const matrix& C_u, const matrix& C_y) {
    int rows = C_du.rows + C_u.rows + C_y.rows;
    int columns = C_du.columns;
    
    matrix A_constrain(rows, columns);
    A_constrain.set_block(0, 0, C_du);
    A_constrain.set_block(C_du.rows, 0, C_u);
    A_constrain.set_block(C_du.rows + C_u.rows, 0, C_y);

    return A_constrain;
}

matrix MPC::construct_constraint_b_du(const matrix& du_max, const matrix& du_min) {
    
    matrix L(h*nu, nu);
    matrix I = mops.eye(nu);
    for (int i = 0;i < h;++i) {
        L.set_block(i*nu,0,I);
    }
    matrix L_du_max = L * du_max;
    matrix L_du_min = L * du_min * -1;
    matrix b_du(2*h*nu, 1);
    b_du.set_block(0, 0, L_du_max);
    b_du.set_block(h*nu, 0, L_du_min);
    return b_du;
}

matrix MPC::construct_constraint_b_u(const matrix& u_max, const matrix& u_min) {
    matrix L(h*nu, nu);
    matrix I = mops.eye(nu);
    for (int i = 0;i < h;++i) {
        L.set_block(i*nu,0,I);
    }
    matrix L_u_max = L * u_max;
    matrix L_u_min = L * u_min * -1;
    matrix b_u(2*h*nu, 1);
    b_u.set_block(0, 0, L_u_max);
    b_u.set_block(h*nu, 0, L_u_min);
    return b_u;
}

matrix MPC::construct_constraint_b_Lu() {
    matrix L(h*nu, nu);
    matrix I = mops.eye(nu);
    for (int i = 0;i < h;++i) {
        L.set_block(i*nu,0,I);
    }
    matrix b_Lu(2*h*nu, nu);
    b_Lu.set_block(0, 0, L);
    b_Lu.set_block(h*nu, 0, L*(-1));
    return b_Lu;
}

matrix MPC::construct_constraint_b_yu() {
    int rows = 2 * P.rows;
    int columns = P.columns;

    matrix b_yu(rows, columns);
    b_yu.set_block(0, 0, P);
    b_yu.set_block(P.rows, 0, P*(-1));
    return b_yu;
}

matrix MPC::construct_constraint_b_y(const matrix& x_max, const matrix& x_min) {
    matrix L(h*nx, nx);
    matrix I = mops.eye(nx);
    for (int i = 0;i < h;++i) {
        L.set_block(i*nx,0,I);
    }
    matrix L_u_max = L * x_max;
    matrix L_u_min = L * x_min * -1;
    matrix b_y(2*h*nx, 1);
    b_y.set_block(0, 0, L_u_max);
    b_y.set_block(h*nx, 0, L_u_min);
    
    return b_y;
}

matrix MPC::construct_b_c(const matrix& b_du, const matrix& b_u, const matrix& b_y) {
    int rows = b_du.rows + b_u.rows + b_y.rows;
    int columns = b_du.columns;
    matrix b_c(rows,columns);
    b_c.set_block(0, 0, b_du);
    b_c.set_block(b_du.rows, 0, b_u);
    b_c.set_block(b_du.rows + b_u.rows, 0, b_y);
    return b_c;
}

matrix MPC::construct_b_cu(const matrix& b_du, const matrix& b_u, const matrix& b_y, const matrix& b_Lu) {
    // b_cu = [O;b_LU;O]
    int rows = b_du.rows + b_u.rows + b_y.rows;
    int columns = b_Lu.columns;
    matrix b_cu(rows,columns);
    matrix O1 = mops.zeros(b_du.rows, columns);
    matrix O3 = mops.zeros(b_y.rows, columns);
    b_cu.set_block(0, 0, O1);
    b_cu.set_block(O1.rows, 0, b_Lu);
    b_cu.set_block(O1.rows + b_Lu.rows, 0, O3);
    return b_cu;
}

matrix MPC::construct_b_cdu(const matrix& b_du, const matrix& b_u, const matrix& b_y, const matrix& b_yu) {
    // b_cdu = [O;O;b_yu]
    int rows = b_du.rows + b_u.rows + b_y.rows;
    int columns = b_yu.columns;
    matrix b_cdu(rows,columns);
    matrix O1 = mops.zeros(b_du.rows, columns);
    matrix O2 = mops.zeros(b_u.rows, columns);
    b_cdu.set_block(0, 0, O1);
    b_cdu.set_block(O1.rows, 0, O2);
    b_cdu.set_block(O1.rows + O2.rows, 0, b_yu);
    return b_cdu;
}

matrix MPC::compute_constraint_b(const matrix& u_prev, const matrix& x_current) {
    // b = b_c - b_cu*u_prev - P * current_z
    // where z is the augmented state z(k) = [x(k);u(k-1)]
    matrix current_z(nx + nu, 1);
    current_z.set_block(0, 0, x_current);
    current_z.set_block(nx, 0 , u_prev);

    int rows = b_c.rows;
    int columns = b_c.columns;
    matrix b(rows, columns);
    b = b_c - b_cu*u_prev - (b_cdu * current_z);
    return b;
}

matrix MPC::construct_c_trans_L1(const matrix& W_x, const matrix& W_u, const matrix& W_du) {
    int nu = W_u.rows; 
    int nx = W_x.rows;
    int columns = h * (nx + 3 * nu);
    matrix c_T(1, columns);
    matrix w_du = mops.diag(W_du).T();
    matrix w_u = mops.diag(W_u).T();
    matrix w_x = mops.diag(W_x).T();

    int offset = 0;
    for (int i = 0; i < h; ++i) {
        c_T.set_block(0, offset, w_du);
        offset += nu;
    }

    for (int i = 0; i < h; ++i) {
        c_T.set_block(0, offset, w_du);
        offset += nu;
    }

    for (int i = 0; i < h; ++i) {
        c_T.set_block(0, offset , w_u);
        offset += nu;
    }

    for (int i = 0; i < h; ++i) {
        c_T.set_block(0, offset, w_x);
        offset += nx;
    }

    return c_T;
}
    
std::vector<ConstraintType> MPC::fill_contraints(const matrix& b) {
    return std::vector<ConstraintType>(b.rows, LEQ);
}

matrix MPC::construct_A_simplex() {
    int A_constrain_rows = A_constraint.rows;
    int u_rows = h * nu;
    int x_rows = h * nx;

    int tot_rows = A_constrain_rows + 2 * u_rows + 2 * x_rows;
    int tot_cols = 3 * u_rows + x_rows; // du+, du-, epsu, epsx

    matrix A_simp = mops.zeros(tot_rows, tot_cols);

    matrix T_u = mops.lower_toeplitz_I(u_rows, nu);
    matrix m_A_constraint = A_constraint * (-1.0);
    matrix m_H = H * (-1.0);
    matrix m_T_u = T_u * (-1.0);

    // Identity matrices for slacks

    matrix m_I_u = mops.eye(u_rows) * (-1.0);
    matrix m_I_x = mops.eye(x_rows) * (-1.0);

    int offset = 0; // offset for rows

    // 1. Physical Constraints [A_orig, -A_orig, 0, 0]
    A_simp.set_block(offset, 0, A_constraint);
    A_simp.set_block(offset, u_rows, m_A_constraint);
    offset += A_constrain_rows;

    // 2. Upper u bound [T_u, -T_u, -I_u, 0]
    A_simp.set_block(offset, 0, T_u);
    A_simp.set_block(offset, u_rows, m_T_u);
    A_simp.set_block(offset, 2 * u_rows, m_I_u);
    offset += u_rows;

    // 3. Lower u bound [-T_u, T_u, -I_u, 0]
    A_simp.set_block(offset, 0, m_T_u);
    A_simp.set_block(offset, u_rows, T_u);
    A_simp.set_block(offset, 2 * u_rows, m_I_u);
    offset += u_rows;

    // 4. Upper x bound (Tracking error) [H, -H, 0, -I_x]
    A_simp.set_block(offset, 0, H);
    A_simp.set_block(offset, u_rows, m_H);
    A_simp.set_block(offset, 3 * u_rows, m_I_x); 
    offset += x_rows;

    // 5. Lower x bound (Tracking error) [-H, H, 0, -I_x]
    A_simp.set_block(offset, 0, m_H);
    A_simp.set_block(offset, u_rows, H);
    A_simp.set_block(offset, 3 * u_rows, m_I_x);

    return A_simp;

}

matrix MPC::construct_b_simplex(const matrix& b_constraint, const matrix& u_past, const matrix& e) {
    // b_simp = [b_constraint;-u_past;u_past;e;-e]
    
    int tot_rows = b_constraint.rows + 2 * h * nu + 2 * h * nx;
    matrix b_simp(tot_rows,1);

    int offset = 0;
    b_simp.set_block(0, 0, b_constraint);
    offset += b_constraint.rows;

    b_simp.set_block(offset, 0, u_past * (-1.0));
    offset += h * nu;

    b_simp.set_block(offset, 0, u_past);
    offset += h * nu;

    b_simp.set_block(offset, 0, e);
    offset += h * nx;

    b_simp.set_block(offset, 0, e * (-1.0));

    return b_simp;
}

matrix MPC::solve(const matrix& measurement, const matrix& reference) {

    current_z.set_block(0, 0, measurement);
    current_z.set_block(nx, 0, prev_control);

    b_constraint = compute_constraint_b(prev_control, measurement);
    b_qp = b_constraint * (-1.0);
    matrix free_resp = P * current_z;

    matrix e = reference - free_resp;

    if (solver_type == LP) {
        for (int i = 0; i < h; ++i) {
        u_past.set_block(i * nu, 0 , prev_control);
        }

        b_simp = construct_b_simplex(b_constraint, u_past, e);

        


        //std::cout << "Initializing simplex" << "\n";
        //std::cout << "Solving" << "\n";
        lp_solver->reset(b_simp);
        lp_solver->solve();
        //std::cout << "returning solution" << "\n";
        matrix Z = lp_solver->return_solution();

        matrix du_pos(nu, 1);
        matrix du_neg(nu, 1);

        //std::cout << "extrackting control" << "\n";
        for (int i = 0; i < nu; ++i) {
            du_pos(i, 0) = Z(i,0);
            du_neg(i, 0) = Z(i + (h * nu),0);
        }

        du = du_pos - du_neg;
    } else if (solver_type == QP) {
        matrix c = H.T() * (-1.0) * W_x_block * e;
        qp_solver->reset(b_qp, c);
        qp_solver->solve();
        
        du = qp_solver->return_solution();
    }


    prev_control = prev_control + du;

    return prev_control;
}

MPC::MPC(const matrix& A, const matrix& B, const matrix& C,
         const matrix& W_x, const matrix& W_u, const matrix& W_du,
         const matrix& du_max, const matrix& du_min, const matrix& u_max,
         const matrix& u_min, const matrix& x_max, const matrix& x_min,
         const matrix& x0, int horizon, SolverType solver_type) {
    
    this->solver_type = solver_type;
    h = horizon;
    A_a = construct_A_a(A, B);  
    B_a = construct_B_a(A, B);
    C_a = construct_C_a(C, B);
    nu = B.columns;
    nx = A.rows;

    P = construct_P(A_a, C_a, horizon);
    //std::cout << "P: " << "\n";
    //mops.print(P);

    P_x = construct_P_x(A_a, horizon);
    //std::cout << "P_x: " << "\n";
    //mops.print(P_x);

    H = construct_H(A_a, B_a, C_a, horizon);
    //std::cout << "H: " << "\n";
    //mops.print(H);

    H_x = construct_H_x(A_a, B_a, horizon);
    //std::cout << "H_x: " << "\n";
    //mops.print(H_x);
    //cost = construct_linear_cost(H_x, W_x, W_u)
    
    C_du = construct_constraint_C_du();
    //std::cout << "C_du: " << "\n";
    //mops.print(C_du);

    C_u = construct_constraint_C_u();
    //std::cout << "C_u: " << "\n";
    //mops.print(C_u);
    C_y = construct_constraint_C_y();
    //std::cout << "C_y: " << "\n";
    //mops.print(C_y);
    
    
    A_constraint = construct_constraint_A(C_du, C_u, C_y);
    //std::cout << "A_constraint: " << "\n";
    //mops.print(A_constraint);
    
    b_du = construct_constraint_b_du(du_max, du_min);
    //std::cout << "b_du: " << "\n";
    //mops.print(b_du);
    std::cout << "test" << "\n";
    b_u = construct_constraint_b_u(u_max, u_min);
    //std::cout << "b_u: " << "\n";
    //mops.print(b_u);

    b_Lu = construct_constraint_b_Lu();
    //std::cout << "b_Lu: " << "\n";
    //mops.print(b_Lu);
    
    b_y = construct_constraint_b_y(x_max, x_min);
    
    b_c =  construct_b_c(b_du, b_u, b_y);
    
    b_cu =  construct_b_cu(b_du, b_u, b_y, b_Lu);

    b_yu = construct_constraint_b_yu();

    b_cdu = construct_b_cdu(b_du, b_u, b_y, b_yu);
    
    matrix u0 = mops.zeros(B.columns,1);
    
    current_z = matrix(nx + nu, 1);
    prev_control = mops.zeros(nu, 1);

    if (solver_type == LP) {
        c_T = construct_c_trans_L1(W_x, W_u, W_du);
        //std::cout << "c_trans: " << "\n";
        //mops.print(c_T);

        A_simp = construct_A_simplex();
        

        
        std::cout << "Amount of A_simp rows: " << A_simp.rows << "\n";
        u_past = matrix(h * nu, 1);
        

        constraint_types = std::vector<ConstraintType>(A_simp.rows, LEQ);

        // Initialize simplex solver 
        matrix dummy_b = mops.zeros(A_simp.rows, 1);
        lp_solver = new Simplex(A_simp, dummy_b, c_T, constraint_types);
    } else if (solver_type == QP) {
        int u_rows = h * nu;
        matrix T_u = mops.lower_toeplitz_I(u_rows, nu);
        W_x_block = mops.create_block_diagonal(W_x, h);
        matrix W_du_block = mops.create_block_diagonal(W_du, h);
        matrix W_u_block = mops.create_block_diagonal(W_u, h);
        matrix Q = H.T() * W_x_block * H + T_u.T() * W_u_block * T_u + W_du_block;
        
        A_qp = A_constraint * (-1.0);
        
        matrix dummy_g = mops.zeros(u_rows, 1);
        matrix dummy_b_constraint = mops.zeros(A_constraint.rows, 1);
        qp_solver = new Active_Set(Q, dummy_g, A_qp, dummy_b_constraint);
    }
    
}


MPC::~MPC() {
    if (solver_type == LP) {
        delete lp_solver;
    } else if (solver_type == QP) {
        delete qp_solver;
    }
}