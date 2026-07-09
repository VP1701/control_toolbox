// test_pid_comparison.cpp
// Runs the same plant/reference scenario as test_mpc_comparison.cpp, but
// controlled with PidVelocity (incremental form) and PidPositional
// (positional form w/ derivative filter + anti-windup) instead of MPC.
// Results are saved to CSV so they can be plotted/compared against
// results_LP.csv / results_QP.csv from the MPC test.

#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <filesystem>
#include "matrix.h"
#include "pid.h"

Matrix mops;

// ---------------------------------------------------------------------
// The controlled output is x1 (index 1 of the state), matching the
// setpoint used in the MPC test (setpoint(1,0)).
// ---------------------------------------------------------------------

void run_simulation_velocity(PidVelocity& controller, const matrix& A, const matrix& B,
                              const matrix& x0, int sim_steps, double u_low, double u_high,
                              const char* filename) {
    controller.reset();

    double r = 10.0;
    matrix x_current = x0;
    std::ofstream outfile(filename);
    outfile << "k,x_r,x0,x1,u\n";

    double tot_ms = 0.0;
    typedef std::chrono::high_resolution_clock Clock;

    for (int k = 0; k < sim_steps; ++k) {
        if (k == 51) {
            r = 7.5;
        }

        double y = x_current(1, 0);
        double e = r - y;

        auto t1 = Clock::now();
        double u_unclamped = controller.calculate_control(e);
        double u = std::clamp(u_unclamped, u_low, u_high);
        auto t2 = Clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
        tot_ms += ns.count() / 1000000.0;

        matrix u_mat(1, 1);
        u_mat(0, 0) = u;
        x_current = (A * x_current) + (B * u_mat);

        outfile << k << ","
                << r << ","
                << x_current(0, 0) << ","
                << x_current(1, 0) << ","
                << u << "\n";
    }

    std::cout << "[PidVelocity] Average time per iteration: "
              << tot_ms / sim_steps << " ms\n";
    outfile.close();
}

void run_simulation_positional(PidPositional& controller, const matrix& A, const matrix& B,
                                const matrix& x0, int sim_steps,
                                const char* filename) {
    controller.reset();

    double r = 10.0;
    matrix x_current = x0;
    std::ofstream outfile(filename);
    outfile << "k,x_r,x0,x1,u\n";

    double tot_ms = 0.0;
    typedef std::chrono::high_resolution_clock Clock;

    for (int k = 0; k < sim_steps; ++k) {
        if (k == 51) {
            r = 7.5;
        }

        double y = x_current(1, 0);

        auto t1 = Clock::now();
        double u = controller.calculate_control(r, y);
        auto t2 = Clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
        tot_ms += ns.count() / 1000000.0;

        matrix u_mat(1, 1);
        u_mat(0, 0) = u;
        x_current = (A * x_current) + (B * u_mat);

        outfile << k << ","
                << r << ","
                << x_current(0, 0) << ","
                << x_current(1, 0) << ","
                << u << "\n";
    }

    std::cout << "[PidPositional] Average time per iteration: "
              << tot_ms / sim_steps << " ms\n";
    outfile.close();
}

int main() {
    int nx = 2;
    int nu = 1;
    int sim_steps = 150;
    double h = 1.0; // sample time, matches the MPC's sample time

    matrix A(nx, nx);
    A(0, 0) = 0.9002;  A(0, 1) = -0.095;
    A(1, 0) = 0.095;   A(1, 1) = 0.9952;

    matrix B(nx, nu);
    B(0, 0) = 0.095;
    B(1, 0) = 0.004833;

    matrix x0 = mops.zeros(nx, 1);

    double u_low = -25.0;
    double u_high = 25.0;

    const std::string out_dir = "control_toolbox/test_results";
    std::filesystem::create_directories(out_dir);

    std::cout << "=== Building PidVelocity controller ===\n";
    PidVelocity vel_controller(
        /*K_p=*/0.5, /*K_d=*/0.0, /*K_i=*/0.05, /*h=*/h,
        ControlMethod::Tustin, /*antiwindup=*/true, u_low, u_high);
    std::cout << "Running PidVelocity simulation...\n";
    run_simulation_velocity(vel_controller, A, B, x0, sim_steps, u_low, u_high,
                             (out_dir + "/results_PID_velocity.csv").c_str());
    std::cout << "PidVelocity done. Results saved to " << out_dir << "/results_PID_velocity.csv\n\n";

    std::cout << "=== Building PidPositional controller ===\n";
    PidPositional pos_controller(
        /*K_p=*/1.0, /*K_i=*/0.05, /*K_d=*/0.05, /*h=*/h,
        u_low, u_high, /*T_f=*/1.0, /*T_aw=*/1.0, /*beta=*/1.0);
    std::cout << "Running PidPositional simulation...\n";
    run_simulation_positional(pos_controller, A, B, x0, sim_steps,
                               (out_dir + "/results_PID_positional.csv").c_str());
    std::cout << "PidPositional done. Results saved to " << out_dir << "/results_PID_positional.csv\n\n";

    std::cout << "Both PID simulations complete. Compare " << out_dir << "/results_PID_velocity.csv, "
              << out_dir << "/results_PID_positional.csv and your MPC results.\n";
    return 0;
}