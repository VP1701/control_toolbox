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
    setpoint(0,0) = 0.0;
    setpoint(1,0) = 10.0;

    matrix reference(horizon * nx, 1);
    for (int i = 0; i < horizon; ++i) {
        reference.set_block(i * nx, 0, setpoint);
    }

    std::ofstream outfile(filename);
    outfile << "k,x_r,x0,x1,u\n";

    matrix x_current = x0;

    double tot_ms = 0.0;
    

    for (int k = 0; k < sim_steps; ++k) {
        if (k == 51) {
            setpoint(0,0) = 0.0;
            setpoint(1,0) = 7.5;
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
        tot_ms += ms;
        x_current = (A * x_current) + (B * u_optimal);

        outfile << k << ","
                << setpoint(1,0) << ","
                << x_current(0,0) << ","
                << x_current(1,0) << ","
                << u_optimal(0,0) << "\n";
    }

    

    std::cout << "Average time taken for 1 MPC iteration " << tot_ms / sim_steps << " ms\n";
    outfile.close();
}

int main() {

    int nx = 2;
    int nu = 1;
    int horizon = 15;
    int horizon_LP = 5;
    int sim_steps = 150;

    matrix A(nx, nx);
    A(0,0) = 0.9002;  A(0,1) = -0.095;
    A(1,0) = 0.095;   A(1,1) = 0.9952;

    matrix B(nx, nu);
    B(0,0) = 0.095;
    B(1,0) = 0.004833;

    matrix C = mops.eye(nx);

    matrix W_x = mops.zeros(nx, nx);
    W_x(0,0) = 0.0;
    W_x(1,1) = 100.0;

    matrix W_x_LP = mops.zeros(nx, nx);
    W_x_LP(0,0) = 0.0;
    W_x_LP(1,1) = 40.0;

    matrix W_u  = mops.eye(nu) * 0.1;
    matrix W_du = mops.eye(nu) * 0.1;

    matrix W_u_LP  = mops.eye(nu) * 0.5;
    matrix W_du_LP = mops.eye(nu) * 1.0;

    matrix du_max = mops.ones(nu, 1) * 10.0;
    matrix du_min = mops.ones(nu, 1) * (-10.0);
    matrix u_max  = mops.ones(nu, 1) * 25.0;
    matrix u_min  = mops.ones(nu, 1) * (-25.0);
    matrix x_max  = mops.ones(nx, 1) * 50.0;
    matrix x_min  = mops.ones(nx, 1) * (-50.0);

    matrix x0 = mops.zeros(nx, 1);
    
    // --- LP test ---
    std::cout << "=== Building LP controller ===\n";
    MPC lp_controller(A, B, C, W_x_LP, W_u_LP, W_du_LP,
                      du_max, du_min, u_max, u_min, x_max, x_min,
                      x0, horizon_LP, LP);
    std::cout << "Running LP simulation...\n";
    run_simulation(lp_controller, A, B, x0, horizon_LP, sim_steps, "results_LP.csv");
    std::cout << "LP done. Results saved to results_LP.csv\n\n";
    
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