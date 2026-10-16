// mpc_quad_hover_test.cpp
// Linear MPC (QP) on a quadrotor translational model linearized at hover.
//
// Model (z-up, psi = 0, Kamel et al. 2017 structure):
//   x = [px, py, pz, vx, vy, vz, phi, theta]      (deviation from hover = absolute)
//   u = [phi_r, theta_r, T - g]                   (deviation from u_bar = [0, 0, g])
//
//   p_dot     = v
//   v_dot     = [g*theta, -g*phi, T - g] - diag(A_drag) * v
//   phi_dot   = (K * phi_r - phi) / tau
//   theta_dot = (K * theta_r - theta) / tau
//
// Parameters: g = 9.81, A_drag = 0.01, K = 1, tau = 2*zeta/wn + Td
//             with zeta = 0.7, wn = 12, Td = 0.015 -> tau = 0.131667 s
// Continuous A, B = Jacobians of mpc_model at x = 0, u = [0, 0, g].
// Discretized with exact ZOH, Ts = 0.1 s:
//   Ad = expm(A*Ts), Bd = int_0^Ts expm(A*s) ds * B
// computed from expm([A B; 0 0]*Ts). Regenerate with
// matlab_testcode/mpc_quad_hover_test.m (c2d(ss(A,B,eye(8),0), 0.1, 'zoh')).
//
// The plant in the simulation is the same linear model -> constraints on
// predicted outputs must hold exactly in closed loop (up to solver tolerance).

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include "matrix.h"
#include "mpc.h"

Matrix mops;

namespace {

constexpr int NX = 8;
constexpr int NU = 3;
constexpr int HORIZON = 20;
constexpr double TS = 0.1;
constexpr double G = 9.81;

constexpr double TILT_MAX = 0.35;        // [rad], ~20 deg, on phi_r/theta_r and phi/theta
constexpr double T_DEV_MIN = -0.7 * G;   // T >= 0.3 g
constexpr double T_DEV_MAX = 0.8 * G;    // T <= 1.8 g
constexpr double CONSTR_TOL = 1.0e-5;    // allowed constraint violation in closed loop

enum StateIdx { PX = 0, PY, PZ, VX, VY, VZ, PHI, THETA };
enum InputIdx { PHI_R = 0, THETA_R, T_DEV };

void build_model(matrix& Ad, matrix& Bd) {
    // Exact ZOH discretization at Ts = 0.1 s, nonzero entries only.
    Ad = mops.zeros(NX, NX);
    Ad(PX, PX) = 1.0;
    Ad(PX, VX) = 9.995001666250e-02;
    Ad(PX, THETA) = 3.865911430106e-02;
    Ad(PY, PY) = 1.0;
    Ad(PY, VY) = 9.995001666250e-02;
    Ad(PY, PHI) = -3.865911430106e-02;
    Ad(PZ, PZ) = 1.0;
    Ad(PZ, VZ) = 9.995001666250e-02;
    Ad(VX, VX) = 9.990004998334e-01;
    Ad(VX, THETA) = 6.868961371220e-01;
    Ad(VY, VY) = 9.990004998334e-01;
    Ad(VY, PHI) = -6.868961371220e-01;
    Ad(VZ, VZ) = 9.990004998334e-01;
    Ad(PHI, PHI) = 4.679032800953e-01;
    Ad(THETA, THETA) = 4.679032800953e-01;

    Bd = mops.zeros(NX, NU);
    Bd(PX, THETA_R) = 1.037453978562e-02;
    Bd(PY, PHI_R) = -1.037453978562e-02;
    Bd(PZ, T_DEV) = 4.998333749917e-03;
    Bd(VX, THETA_R) = 2.936135263372e-01;
    Bd(VY, PHI_R) = -2.936135263372e-01;
    Bd(VZ, T_DEV) = 9.995001666250e-02;
    Bd(PHI, PHI_R) = 5.320967199047e-01;
    Bd(THETA, THETA_R) = 5.320967199047e-01;
}

struct Metrics {
    double max_abs_tilt = 0.0;       // max |phi|, |theta| (states)
    double max_abs_tilt_ref = 0.0;   // max |phi_r|, |theta_r| (inputs)
    double min_t_dev = 0.0;
    double max_t_dev = 0.0;
    double max_abs_u = 0.0;
    double first_theta_r = 0.0;      // sign check for lateral step
    double avg_solve_ms = 0.0;
    double max_solve_ms = 0.0;
    double avg_iters = 0.0;
    int max_iters = 0;
    double ms_per_iter_worst = 0.0;  // solve time / iterations at the max-iteration step
    matrix x_final;
};

void fill_reference(matrix& reference, const matrix& setpoint) {
    // constant setpoint over the whole horizon
    for (int i = 0; i < HORIZON; ++i) {
        reference.set_block(i * NX, 0, setpoint);
    }
}

// setpoint_0 is used for k < change_k, setpoint_1 from k = change_k on.
// Pass change_k >= sim_steps for a constant setpoint.
Metrics run_scenario(MPC& controller, const matrix& Ad, const matrix& Bd,
                     const matrix& setpoint_0, const matrix& setpoint_1,
                     int change_k, int sim_steps, const char* filename) {
    matrix setpoint = setpoint_0;
    matrix reference(HORIZON * NX, 1);
    fill_reference(reference, setpoint);

    std::ofstream outfile(filename);
    outfile << "k,t,px,py,pz,vx,vy,vz,phi,theta,phi_r,theta_r,T,"
               "px_ref,py_ref,pz_ref,solve_ms,iters\n";

    Metrics m;
    matrix x = mops.zeros(NX, 1);
    double total_ms = 0.0;
    int total_iters = 0;

    for (int k = 0; k < sim_steps; ++k) {
        if (k == change_k) {
            setpoint = setpoint_1;
            fill_reference(reference, setpoint);
        }
        const auto t_start = std::chrono::steady_clock::now();
        matrix u = controller.solve(x, reference);
        const auto t_end = std::chrono::steady_clock::now();
        const double solve_ms =
            std::chrono::duration<double, std::milli>(t_end - t_start).count();
        total_ms += solve_ms;
        m.max_solve_ms = std::fmax(m.max_solve_ms, solve_ms);

        const int iters = controller.return_iters();
        total_iters += iters;
        if (iters > m.max_iters) {
            m.max_iters = iters;
            m.ms_per_iter_worst = solve_ms / iters;
        }

        x = (Ad * x) + (Bd * u);

        if (k == 0) {
            m.first_theta_r = u(THETA_R, 0);
        }
        m.max_abs_tilt = std::fmax(m.max_abs_tilt,
                                   std::fmax(std::fabs(x(PHI, 0)), std::fabs(x(THETA, 0))));
        m.max_abs_tilt_ref = std::fmax(m.max_abs_tilt_ref,
                                       std::fmax(std::fabs(u(PHI_R, 0)), std::fabs(u(THETA_R, 0))));
        m.min_t_dev = std::fmin(m.min_t_dev, u(T_DEV, 0));
        m.max_t_dev = std::fmax(m.max_t_dev, u(T_DEV, 0));
        for (int i = 0; i < NU; ++i) {
            m.max_abs_u = std::fmax(m.max_abs_u, std::fabs(u(i, 0)));
        }

        outfile << k << "," << (k + 1) * TS;
        for (int i = 0; i < NX; ++i) {
            outfile << "," << x(i, 0);
        }
        outfile << "," << u(PHI_R, 0) << "," << u(THETA_R, 0) << "," << u(T_DEV, 0) + G
                << "," << setpoint(PX, 0) << "," << setpoint(PY, 0) << "," << setpoint(PZ, 0)
                << "," << solve_ms << "," << iters << "\n";
    }

    m.avg_solve_ms = total_ms / sim_steps;
    m.avg_iters = static_cast<double>(total_iters) / sim_steps;
    m.x_final = x;
    return m;
}

int n_failed = 0;

void print_timing(const Metrics& m) {
    std::cout << "  solve time: avg " << m.avg_solve_ms << " ms, worst "
              << m.max_solve_ms << " ms\n";
    std::cout << "  iterations: avg " << m.avg_iters << ", worst "
              << m.max_iters << "\n";
    if (m.max_iters > 0) {
        std::cout << "  time per iteration (worst step): "
                  << m.ms_per_iter_worst << " ms\n";
    }
}

void check(bool condition, const char* description, double value) {
    std::cout << (condition ? "  PASS  " : "  FAIL  ") << description
              << "  (value: " << value << ")\n";
    if (!condition) {
        ++n_failed;
    }
}

}  // namespace

int main() {
    matrix Ad;
    matrix Bd;
    build_model(Ad, Bd);
    matrix C = mops.eye(NX);

    // Weights (velocity-form MPC: W_du penalizes input rate, W_u set to 0,
    // see note on W_u in the accompanying review)
    matrix W_x = mops.zeros(NX, NX);
    W_x(PX, PX) = 10.0;
    W_x(PY, PY) = 10.0;
    W_x(PZ, PZ) = 10.0;
    W_x(VX, VX) = 1.0;
    W_x(VY, VY) = 1.0;
    W_x(VZ, VZ) = 1.0;
    W_x(PHI, PHI) = 1.0;
    W_x(THETA, THETA) = 1.0;

    matrix W_u = mops.zeros(NU, NU);

    matrix W_du = mops.zeros(NU, NU);
    W_du(PHI_R, PHI_R) = 5.0;
    W_du(THETA_R, THETA_R) = 5.0;
    W_du(T_DEV, T_DEV) = 0.05;

    // Input rate limits per step (Ts = 0.1 s)
    matrix du_max(NU, 1);
    du_max(PHI_R, 0) = 0.2;
    du_max(THETA_R, 0) = 0.2;
    du_max(T_DEV, 0) = 5.0;
    matrix du_min = du_max * (-1.0);

    // Input limits (deviation from hover)
    matrix u_max(NU, 1);
    u_max(PHI_R, 0) = TILT_MAX;
    u_max(THETA_R, 0) = TILT_MAX;
    u_max(T_DEV, 0) = T_DEV_MAX;
    matrix u_min(NU, 1);
    u_min(PHI_R, 0) = -TILT_MAX;
    u_min(THETA_R, 0) = -TILT_MAX;
    u_min(T_DEV, 0) = T_DEV_MIN;

    // Output limits: wide on position and velocity, tilt limit on phi/theta states
    matrix x_max = mops.ones(NX, 1) * 100.0;
    x_max(VX, 0) = 10.0;
    x_max(VY, 0) = 10.0;
    x_max(VZ, 0) = 10.0;
    x_max(PHI, 0) = TILT_MAX;
    x_max(THETA, 0) = TILT_MAX;
    matrix x_min = x_max * (-1.0);

    matrix x0 = mops.zeros(NX, 1);
    const int sim_steps = 100;  // 10 s

    // --- Scenario 1: hover hold --------------------------------------------
    std::cout << "=== Scenario 1: hover hold (r = 0) ===\n";
    {
        MPC mpc(Ad, Bd, C, W_x, W_u, W_du, du_max, du_min, u_max, u_min,
                x_max, x_min, x0, HORIZON, QP);
        matrix sp = mops.zeros(NX, 1);
        Metrics m = run_scenario(mpc, Ad, Bd, sp, sp, sim_steps, sim_steps,
                                 "results_quad_hover.csv");
        print_timing(m);
        check(m.max_abs_u < 1.0e-9, "u stays at hover (max |u_dev|)", m.max_abs_u);
        check(m.x_final.L2() < 1.0e-9, "state stays at origin (||x||)", m.x_final.L2());
    }

    // --- Scenario 2: altitude step pz = 2 m ----------------------------------
    std::cout << "=== Scenario 2: altitude step (pz_ref = 2 m) ===\n";
    {
        MPC mpc(Ad, Bd, C, W_x, W_u, W_du, du_max, du_min, u_max, u_min,
                x_max, x_min, x0, HORIZON, QP);
        matrix sp = mops.zeros(NX, 1);
        sp(PZ, 0) = 2.0;
        Metrics m = run_scenario(mpc, Ad, Bd, sp, sp, sim_steps, sim_steps,
                                 "results_quad_alt.csv");
        print_timing(m);
        check(std::fabs(m.x_final(PZ, 0) - 2.0) < 0.02, "pz settles to 2 m (|error|)",
              std::fabs(m.x_final(PZ, 0) - 2.0));
        check(m.max_abs_tilt < 1.0e-9, "no tilt excited (max |phi|,|theta|)", m.max_abs_tilt);
        check(m.min_t_dev >= T_DEV_MIN - CONSTR_TOL, "T >= 0.3 g (min T - g)", m.min_t_dev);
        check(m.max_t_dev <= T_DEV_MAX + CONSTR_TOL, "T <= 1.8 g (max T - g)", m.max_t_dev);
    }

    // --- Scenario 3: lateral step px = 5 m, tilt limit active ----------------
    std::cout << "=== Scenario 3: lateral step (px_ref = 5 m) ===\n";
    {
        MPC mpc(Ad, Bd, C, W_x, W_u, W_du, du_max, du_min, u_max, u_min,
                x_max, x_min, x0, HORIZON, QP);
        matrix sp = mops.zeros(NX, 1);
        sp(PX, 0) = 5.0;
        Metrics m = run_scenario(mpc, Ad, Bd, sp, sp, sim_steps, sim_steps,
                                 "results_quad_lat.csv");
        print_timing(m);
        check(m.first_theta_r > 0.0, "first move pitches toward +x (theta_r > 0, z-up)",
              m.first_theta_r);
        check(std::fabs(m.x_final(PX, 0) - 5.0) < 0.05, "px settles to 5 m (|error|)",
              std::fabs(m.x_final(PX, 0) - 5.0));
        check(std::fabs(m.x_final(PY, 0)) < 1.0e-9, "no motion in y (|py|)",
              std::fabs(m.x_final(PY, 0)));
        check(std::fabs(m.x_final(PZ, 0)) < 1.0e-9, "no motion in z (|pz|)",
              std::fabs(m.x_final(PZ, 0)));
        check(m.max_abs_tilt <= TILT_MAX + CONSTR_TOL, "tilt state within limit (max)",
              m.max_abs_tilt);
        check(m.max_abs_tilt_ref <= TILT_MAX + CONSTR_TOL, "tilt ref within limit (max)",
              m.max_abs_tilt_ref);
    }

    // --- Scenario 4: hover, then px step 5 m at k = 50 ------------------------
    // Same maneuver as scenario 3, but the step arrives mid-run when the
    // warm start is settled at hover -> worst case for a reference change.
    std::cout << "=== Scenario 4: mid-run lateral step (px_ref = 5 m at t = 5 s) ===\n";
    {
        MPC mpc(Ad, Bd, C, W_x, W_u, W_du, du_max, du_min, u_max, u_min,
                x_max, x_min, x0, HORIZON, QP);
        matrix sp_0 = mops.zeros(NX, 1);
        matrix sp_1 = mops.zeros(NX, 1);
        sp_1(PX, 0) = 5.0;
        const int change_k = 50;
        const int steps = 150;  // 15 s
        Metrics m = run_scenario(mpc, Ad, Bd, sp_0, sp_1, change_k, steps,
                                 "results_quad_lat_mid.csv");
        print_timing(m);
        check(std::fabs(m.x_final(PX, 0) - 5.0) < 0.05, "px settles to 5 m (|error|)",
              std::fabs(m.x_final(PX, 0) - 5.0));
        check(m.max_abs_tilt <= TILT_MAX + CONSTR_TOL, "tilt state within limit (max)",
              m.max_abs_tilt);
        check(m.max_abs_tilt_ref <= TILT_MAX + CONSTR_TOL, "tilt ref within limit (max)",
              m.max_abs_tilt_ref);
    }

    std::cout << "\n" << (n_failed == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << n_failed << " failed)\n";
    return (n_failed == 0) ? 0 : 1;
}