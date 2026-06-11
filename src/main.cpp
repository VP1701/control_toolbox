// test_mpc_comparison.cpp
// Runs the same MPC problem with LP and QP solvers and saves results to CSV

#include <iostream>
#include <fstream>
#include <vector>
#include "matrix.h"
#include "simplex.h"
#include "mpc.h"

Matrix mops;

void run_simulation(MPC& controller, const matrix& A, const matrix& B,
                    const matrix& x0, int horizon, int sim_steps,
                    const char* filename) {

    int nx = A.rows;
    int nu = B.columns;

    matrix setpoint(nx, 1);
    setpoint(0,0) = 8.0;
    setpoint(1,0) = 6.0;
    setpoint(2,0) = 0.0;
    setpoint(3,0) = 0.0;

    matrix reference(horizon * nx, 1);
    for (int i = 0; i < horizon; ++i) {
        reference.set_block(i * nx, 0, setpoint);
    }

    std::ofstream outfile(filename);
    outfile << "k,"
            << "r1,r2,"
            << "x1,x2,x3,x4,"
            << "u1,u2\n";

    matrix x_current = x0;

    double tot_ms = 0.0;
    

    for (int k = 0; k < sim_steps; ++k) {
        if (k == 100) {
            setpoint(0,0) = 4.0;
            setpoint(1,0) = 10.0;
            for (int i = 0; i < horizon; ++i) {
                reference.set_block(i * nx, 0, setpoint);
            }
        }
        typedef std::chrono::high_resolution_clock Clock;
        auto t1 = Clock::now();
        matrix u_optimal = controller.solve(x_current, reference);
        auto t2 = Clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1);
        double ms = ns.count() / 1000000.0;
        std::cout << "u_optimal:\n";
        mops.print(u_optimal);
        tot_ms += ms;
        x_current = (A * x_current) + (B * u_optimal);

        outfile << k << ","
                << setpoint(0,0) << ","
                << setpoint(1,0) << ","
                << x_current(0,0) << ","
                << x_current(1,0) << ","
                << x_current(2,0) << ","
                << x_current(3,0) << ","
                << u_optimal(0,0) << ","
                << u_optimal(1,0)
                << "\n";
    }

    

    std::cout << "Average time taken for 1 MPC iteration " << tot_ms / sim_steps << " ms\n";
    outfile.close();
}

int main() {

    int nx = 4;
    int nu = 2;
    int horizon = 20;
    int sim_steps = 300;

    matrix A(nx, nx);
    A(0,0) = 0.97;  A(0,1) = 0.00;  A(0,2) = 0.03;  A(0,3) = 0.00;
    A(1,0) = 0.00;  A(1,1) = 0.96;  A(1,2) = 0.00;  A(1,3) = 0.04;
    A(2,0) = 0.00;  A(2,1) = 0.00;  A(2,2) = 0.95;  A(2,3) = 0.00;
    A(3,0) = 0.00;  A(3,1) = 0.00;  A(3,2) = 0.00;  A(3,3) = 0.94;

    matrix B(nx, nu);
    B(0,0) = 0.08;  B(0,1) = 0.02;
    B(1,0) = 0.01;  B(1,1) = 0.09;
    B(2,0) = 0.05;  B(2,1) = 0.00;
    B(3,0) = 0.00;  B(3,1) = 0.04;

    matrix C = mops.eye(nx);

    matrix W_x = mops.zeros(nx, nx);
    W_x(0,0) = 100.0;
    W_x(1,1) = 100.0;
    W_x(2,2) = 10.0;
    W_x(3,3) = 10.0;


    matrix W_u  = mops.eye(nu) * 0.1;
    matrix W_du = mops.eye(nu) * 1.0;

    matrix u_max(2,1);
    u_max(0,0) = 10.0;
    u_max(1,0) = 10.0;

    matrix u_min(2,1);
    u_min(0,0) = -10.0;
    u_min(1,0) = -10.0;

    matrix du_max(2,1);
    du_max(0,0) = 10.0;
    du_max(1,0) = 10.0;

    matrix du_min(2,1);
    du_min(0,0) = -1.0;
    du_min(1,0) = -1.0;

    matrix x_max(4,1);
    x_max(0,0)=20.0;
    x_max(1,0)=20.0;
    x_max(2,0)=20.0;
    x_max(3,0)=20.0;

    matrix x_min(4,1);
    x_min(0,0)=0.0;
    x_min(1,0)=0.0;
    x_min(2,0)=0.0;
    x_min(3,0)=0.0;

    matrix x0(nx, 1);
    x0(0,0)=1.0;
    x0(1,0)=2.0;
    x0(2,0)=3.0;
    x0(3,0)=4.0;

    
    // --- QP test ---
    std::cout << "=== Building QP controller ===\n";
    MPC qp_controller(A, B, C, W_x, W_u, W_du,
                      du_max, du_min, u_max, u_min, x_max, x_min,
                      x0, horizon, QP);
    std::cout << "Running QP simulation...\n";
    run_simulation(qp_controller, A, B, x0, horizon, sim_steps, "results_QP.csv");
    std::cout << "QP done. Results saved to results_QP.csv\n\n";

    std::cout << "Both simulations complete. Compare results_LP.csv and results_QP.csv.\n";

    return 0;
}