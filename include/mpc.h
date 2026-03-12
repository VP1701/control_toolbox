// mpc.h

#ifndef MPC_H
#define MPC_H

#include "matrix.h"
#include "simplex.h"
#include <iostream>
#include <vector> // needed for constraint type vector

class MPC {
    /* MPC solver
    predictions based on state space models
    state space models are axpected to be stable, 
    no quarantees that MPC stabilizes systems.
    
    */

    private:
        mutable Matrix mops;

        matrix x0; // Initial conditions
        matrix u0; // initial controls. zeros

        int r; // rows (states)
        int c; // columns ()
        int h; // horizon lenght
        int nu; // amount of control input
        int nx; // amount of states

        matrix du_max;
        matrix du_min;
        matrix u_max;
        matrix u_min;
        matrix x_max;
        matrix x_min;

        // Augmented state space model matrices for incremental control
        matrix A_a;
        matrix B_a;
        matrix C_a;

        matrix C_du;
        matrix C_u;
        matrix C_y;

        matrix b_du;
        matrix b_u;
        matrix b_y;
        matrix b_Lu;
        matrix b_yu;
        matrix b_yy;
        
        matrix b_c;
        matrix b_cu;
        matrix b_cdu;

        matrix A_system;
        matrix B_system;

        matrix b;


        // prediction matrices
        matrix P; // state
        matrix P_x; // control of previous time step
        matrix H; // control
        matrix H_x; // control increment
        matrix c_T;

        matrix W_x; // cost for error/state deviation
        matrix W_u; // cost for control magnitude
        matrix W_du; // cost for control change

        matrix A_constraint;
        matrix b_constraint;

        matrix A_simp;
        matrix b_simp;

        // holder for constraints. needed for simplex
        std::vector<ConstraintType> constraint_types;
        // function to fill the constraint type vector
        std::vector<ConstraintType> fill_contraints(const matrix& b);

        // change these to usethe matrices internelly so no need for return or input matrices

        matrix construct_A_a(const matrix& A, const matrix& B);
        matrix construct_B_a(const matrix& A, const matrix& B);
        matrix construct_C_a(const matrix& C, const matrix& B);

        matrix construct_P(const matrix& A_sys, const matrix& C_sys, int horizon);
        matrix construct_P_x(const matrix& A_sys, int horizon);
        
        matrix construct_H(const matrix& A_sys,const matrix& B_sys, const matrix& C_sys, int horizon);
        matrix construct_H_x(const matrix& A_sys,const matrix& B_sys, int horizon);

        // make construct functions for blockdiagonal W_x and W_u

        matrix construct_linear_cost_c(const matrix& H_x, const matrix& W_x, const matrix& W_u);
        matrix construct_constraint_A(const matrix& C_du, const matrix& C_u, const matrix& C_y);
        matrix construct_constraint_b(const matrix& b_du, const matrix& b_u, const matrix& b_y);

        matrix construct_constraint_C_du();
        matrix construct_constraint_C_u();
        matrix construct_constraint_C_y();

        // helpers for helpers of b calculation
        matrix construct_constraint_b_du(const matrix& du_max, const matrix& du_min);
        matrix construct_constraint_b_u(const matrix& u_max, const matrix& u_min);
        matrix construct_constraint_b_Lu();
        matrix construct_constraint_b_yu();
        matrix construct_constraint_b_y(const matrix& x_max, const matrix& x_min);
        //matrix construct_constraint_b_yy(); only needed for CARIMA
        


        // helpers for b calculation
        matrix construct_b_c(const matrix& b_du, const matrix& b_u, const matrix& b_y);
        matrix construct_b_cu(const matrix& b_du, const matrix& b_u, const matrix& b_y, const matrix& b_Lu);
        matrix construct_b_cdu(const matrix& b_du, const matrix& b_u, const matrix& b_y, const matrix& b_yu);

        matrix construct_A_simplex();
        matrix construct_b_simplex();

        // Used to compute b online
        matrix compute_constraint_b(const matrix& u_prev, const matrix& x_current);
        matrix construct_c_trans_L1(const matrix& W_x, const matrix& W_u, const matrix& W_du);


    public:
        MPC(const matrix& A, const matrix& B, const matrix& C,
            const matrix& W_x, const matrix& W_u, const matrix& W_du,
            const matrix& du_max, const matrix& du_min, const matrix& u_max,
            const matrix& u_min, const matrix& x_max, const matrix& x_min,
            const matrix& x0, int horizon);

        // last thing to implement
        matrix solve(const matrix& measurement, const matrix& reference);
};

#endif // MPC_H