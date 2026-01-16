// mpc.cpp
#include "mpc.h"


matrix MPC::construct_S_x(const matrix& A_sys, int horizon) {
    int rows = A_sys.rows;
    int columns = A_sys.columns;
    matrix S_x(horizon * rows, columns);
    matrix current_A = mops.eye(r);
    S_x.set_block(0, 0, current_A);
    
    for (int i = 1;i < horizon; ++i) {
        current_A = current_A*A_sys;
        S_x.set_block(i * rows, 0, current_A);
    }

    return S_x;
}

matrix MPC::construct_S_u(const matrix& A_sys,const matrix& B_sys, int horizon) {
    int rows = horizon * A_sys.rows;
    int columns = horizon * B_sys.columns;
    matrix S_u(rows, columns);

    for (int i = 2; i < rows; ++i) {
        for (int j = 0; j< columns; ++j) {
            S_u.set_block()
        }
    }
}