// nmpc.h

#ifndef NMPC_H
#define NMPC_H

#include "matrix.h"
#include "simplex.h"
#include "active_set.h"
#include <iostream>
#include <vector>

class NMPC {

    private:
        Active_Set* qp_solver;

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

        matrix W_x; // cost for error/state deviation
        matrix W_u; // cost for control magnitude
        matrix W_du; // cost for control change

        matrix prev_control;
        matrix current_z;
        matrix u_past;
        matrix du;
    public:
        // Todo
        NMPC();

        matrix solve(const matrix& measurement, const matrix& reference);
        ~NMPC();
};

#endif // NMPC_H