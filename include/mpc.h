// mpc.h

#ifndef MPC_H
#define MPC_H

#include "matrix.h"
#include "simplex.h"
#include <iostream>


class MPC {
    private:
        mutable Matrix mops;

        matrix x0; // Initial conditions
        
        int r; // rows (states)
        int c; // columns ()
        int h; // horizon lenght

        matrix A_system;
        matrix B_system;
        // prediction matrices
        matrix S_x;
        matrix S_u;

        matrix construct_S_x(const matrix& A_sys, int horizon);
        matrix construct_S_u(const matrix& A_sys,const matrix& B_sys, int horizon);


    public:
        ConstructMPCProblem(const matrix& A, const matrix& B, int horizon);

        matrix solve(const matrix& measurement, const matrix& reference);
}


#endif // MPC_H