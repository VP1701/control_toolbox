// linear_drone.cpp

/*
state vector x = [px, py, pz, vx, vy, vz, phi, theta, T] ^ T




*/
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include "matrix.h"
#include "mpc.h"
#include <array>
#include <cstddef>
#include <limits>

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


void build_model(matrix& Ad, matrix& Bd) {
    // ZOH discretization at Ts = 0.1 s, nonzero entries only. values not in full precision
    Ad = mops.zeros(NX, NX);
    Ad(0,0) = 1.0;
    Ad(1,1) = 1.0;
    Ad(2,2) = 1.0;
    Ad(3,3) = 0.999;
    Ad(4,4) = 0.999;
    Ad(5,5) = 0.999;
    Ad(6,6) = 0.4541;
    Ad(7,7) = 0.4541;

    Ad(0,3) = 0.09995;
    Ad(1,4) = 0.09995;
    Ad(2,5) = 0.09995;
    Ad(4,6) = -0.678;
    Ad(3,7) = 0.678;

    Ad(0,7) = 0.03832;
    Ad(1,6) = -0.03832;

    Bd = mops.zeros(NX, NU);
    Bd(0,1) = 0.01071;
    Bd(1,0) = -0.01071;
    Bd(2,2) = 0.004998;

    Bd(3,1) = 0.3025;
    Bd(4,0) = -0.3025;
    Bd(5,2) = 0.09995;

    Bd(6,0) = 0.5459;
    Bd(7,1) = 0.5459;
    
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

void step_and_log(MPC& controller, const matrix& Ad, const matrix& Bd,
                  const matrix& reference, const matrix& setpoint, int k,
                  matrix& x, Metrics& m, double& total_ms, int& total_iters,
                  std::ofstream& outfile) {
    const auto t_start = std::chrono::steady_clock::now();
    const matrix u = controller.solve(x, reference);
    const auto t_end = std::chrono::steady_clock::now();
    const double solve_ms =
        std::chrono::duration<double, std::milli>(t_end - t_start).count();
    total_ms += solve_ms;
    m.max_solve_ms = std::fmax(m.max_solve_ms, solve_ms);

    const int iters = controller.return_iters();
    total_iters += iters;
    if (iters > m.max_iters) {
        m.max_iters = iters;
        m.ms_per_iter_worst = solve_ms / static_cast<double>(iters);
    }

    // plant update 
    x = (Ad * x) + (Bd * u);

    if (k == 0) {
        m.first_theta_r = u(1, 0);
    }
    m.max_abs_tilt = std::fmax(m.max_abs_tilt,
                               std::fmax(std::fabs(x(6, 0)), std::fabs(x(7, 0))));
    m.max_abs_tilt_ref = std::fmax(m.max_abs_tilt_ref,
                                   std::fmax(std::fabs(u(0, 0)), std::fabs(u(1, 0))));
    m.min_t_dev = std::fmin(m.min_t_dev, u(2, 0));
    m.max_t_dev = std::fmax(m.max_t_dev, u(2, 0));
    for (int i = 0; i < NU; ++i) {
        m.max_abs_u = std::fmax(m.max_abs_u, std::fabs(u(i, 0)));
    }

    // logged state is x(k+1), so the time stamp is (k+1)*Ts
    outfile << k << "," << (k + 1) * TS;
    for (int i = 0; i < NX; ++i) {
        outfile << "," << x(i, 0);
    }
    outfile << "," << u(0, 0) << "," << u(1, 0) << "," << u(2, 0) + G
            << "," << setpoint(0, 0) << "," << setpoint(1, 0) << "," << setpoint(2, 0)
            << "," << solve_ms << "," << iters << "\n";
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
        step_and_log(controller, Ad, Bd, reference, setpoint, k,
                 x, m, total_ms, total_iters, outfile);
    }

    m.avg_solve_ms = total_ms / sim_steps;
    m.avg_iters = static_cast<double>(total_iters) / sim_steps;
    m.x_final = x;
    return m;
}


// --- time-parameterized trajectory scenario ---------------------------------
constexpr std::size_t N_PTS = 6;
using PathPoints = std::array<std::array<double, 3>, N_PTS>;

// Straight segments between the points, flown at constant speed.
struct Trajectory {
    PathPoints points{};
    std::array<double, N_PTS> t_arrive{};  // arrival time at each point [s]
};

Trajectory make_trajectory(const PathPoints& points, double speed) {
    // calculates the trajectory for given setpoints and a constant speed
    Trajectory traj;
    traj.points = points;
    traj.t_arrive[0] = 0.0;
    for (std::size_t j = 0; j + 1U < N_PTS; ++j) {
        const double length = std::hypot(points[j + 1U][0] - points[j][0],
                                          points[j + 1U][1] - points[j][1],
                                          points[j + 1U][2] - points[j][2]);
        traj.t_arrive[j + 1U] = traj.t_arrive[j] + (length / speed);
    }
    return traj;
}

// NX x 1 setpoint on the trajectory at time t: position and velocity
// reference, angles zero. Before t = 0 the start point is held, after the
// last arrival time the end point is held with zero velocity.
matrix trajectory_point(const Trajectory& traj, double t) {
    matrix sp = mops.zeros(NX, 1);
    const std::size_t last = N_PTS - 1U;
    if (t <= 0.0) {
        for (std::size_t a = 0; a < 3U; ++a) {
            sp(static_cast<int>(a), 0) = traj.points[0][a];
        }
        return sp;
    }
    if (t >= traj.t_arrive[last]) {
        for (std::size_t a = 0; a < 3U; ++a) {
            sp(static_cast<int>(a), 0) = traj.points[last][a];
        }
        return sp;
    }

    // segment j: t_arrive[j] <= t < t_arrive[j + 1]
    std::size_t j = 0;
    while ((j + 1U < last) && (t >= traj.t_arrive[j + 1U])) {
        ++j;
    }
    const double seg_time = traj.t_arrive[j + 1U] - traj.t_arrive[j];
    const double s = (t - traj.t_arrive[j]) / seg_time;  // 0..1 along segment
    for (std::size_t a = 0; a < 3U; ++a) {
        const double delta = traj.points[j + 1U][a] - traj.points[j][a];
        const int row = static_cast<int>(a);
        sp(row, 0) = traj.points[j][a] + (s * delta);   // position
        sp(row + 3, 0) = delta / seg_time;              // velocity
    }
    return sp;
}

struct TrajectoryMetrics {
    Metrics base;
    double max_track_err = 0.0;                 // max ||p(t) - p_ref(t)|| [m]
    std::array<double, N_PTS> min_gate_dist{};  // closest approach per point [m]
};

// Tracks the trajectory. Reference block i holds the trajectory at
// t = (k + i + 1) * Ts, matching the MPC prediction x(k+1) ... x(k+N).
TrajectoryMetrics run_trajectory(MPC& controller, const matrix& Ad, const matrix& Bd,
                                 const Trajectory& traj, int sim_steps,
                                 const char* filename) {
    matrix reference(HORIZON * NX, 1);

    std::ofstream outfile(filename);
    outfile << "k,t,px,py,pz,vx,vy,vz,phi,theta,phi_r,theta_r,T,"
               "px_ref,py_ref,pz_ref,solve_ms,iters\n";

    TrajectoryMetrics tm;
    tm.min_gate_dist.fill(std::numeric_limits<double>::max());
    matrix x = mops.zeros(NX, 1);
    double total_ms = 0.0;
    int total_iters = 0;

    for (int k = 0; k < sim_steps; ++k) {
        for (int i = 0; i < HORIZON; ++i) {
            const double t_i = static_cast<double>(k + i + 1) * TS;
            reference.set_block(i * NX, 0, trajectory_point(traj, t_i));
        }
        // reference at t = (k+1) * Ts, same time as the logged state x(k+1)
        const matrix setpoint = trajectory_point(traj, static_cast<double>(k + 1) * TS);

        step_and_log(controller, Ad, Bd, reference, setpoint, k,
                     x, tm.base, total_ms, total_iters, outfile);

        const double track_err = std::hypot(x(0, 0) - setpoint(0, 0),
                                            x(1, 0) - setpoint(1, 0),
                                            x(2, 0) - setpoint(2, 0));
        tm.max_track_err = std::fmax(tm.max_track_err, track_err);
        for (std::size_t j = 0; j < N_PTS; ++j) {
            const double dist = std::hypot(x(0, 0) - traj.points[j][0],
                                           x(1, 0) - traj.points[j][1],
                                           x(2, 0) - traj.points[j][2]);
            tm.min_gate_dist[j] = std::fmin(tm.min_gate_dist[j], dist);
        }
    }

    tm.base.avg_solve_ms = total_ms / sim_steps;
    tm.base.avg_iters = static_cast<double>(total_iters) / sim_steps;
    tm.base.x_final = x;
    return tm;
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
    W_x(0, 0) = 10.0;
    W_x(1, 1) = 10.0;
    W_x(2, 2) = 10.0;
    W_x(3, 3) = 1.0;
    W_x(4, 4) = 1.0;
    W_x(5, 5) = 1.0;
    W_x(6, 6) = 1.0;
    W_x(7, 7) = 1.0;

    matrix W_u = mops.zeros(NU, NU);

    matrix W_du = mops.zeros(NU, NU);
    W_du(0, 0) = 5.0;
    W_du(1, 1) = 5.0;
    W_du(2, 2) = 0.05;

    // Input rate limits per step (Ts = 0.1 s)
    matrix du_max(NU, 1);
    du_max(0, 0) = 0.2;
    du_max(1, 0) = 0.2;
    du_max(2, 0) = 5.0;
    matrix du_min = du_max * (-1.0);

    // Input limits (deviation from hover)
    matrix u_max(NU, 1);
    u_max(0, 0) = TILT_MAX;
    u_max(1, 0) = TILT_MAX;
    u_max(2, 0) = T_DEV_MAX;
    matrix u_min(NU, 1);
    u_min(0, 0) = -TILT_MAX;
    u_min(1, 0) = -TILT_MAX;
    u_min(2, 0) = T_DEV_MIN;

    // Output limits: wide on position and velocity, tilt limit on phi/theta states
    matrix x_max = mops.ones(NX, 1) * 100.0;
    x_max(3, 0) = 10.0;
    x_max(4, 0) = 10.0;
    x_max(5, 0) = 10.0;
    x_max(6, 0) = TILT_MAX;
    x_max(7, 0) = TILT_MAX;
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
        sp(2, 0) = 2.0;
        Metrics m = run_scenario(mpc, Ad, Bd, sp, sp, sim_steps, sim_steps,
                                 "results_quad_alt.csv");
        print_timing(m);
        check(std::fabs(m.x_final(2, 0) - 2.0) < 0.02, "pz settles to 2 m (|error|)",
              std::fabs(m.x_final(2, 0) - 2.0));
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
        sp(0, 0) = 5.0;
        Metrics m = run_scenario(mpc, Ad, Bd, sp, sp, sim_steps, sim_steps,
                                 "results_quad_lat.csv");
        print_timing(m);
        check(m.first_theta_r > 0.0, "first move pitches toward +x (theta_r > 0, z-up)",
              m.first_theta_r);
        check(std::fabs(m.x_final(0, 0) - 5.0) < 0.05, "px settles to 5 m (|error|)",
              std::fabs(m.x_final(0, 0) - 5.0));
        check(std::fabs(m.x_final(1, 0)) < 1.0e-9, "no motion in y (|py|)",
              std::fabs(m.x_final(1, 0)));
        check(std::fabs(m.x_final(2, 0)) < 1.0e-9, "no motion in z (|pz|)",
              std::fabs(m.x_final(2, 0)));
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
        sp_1(0, 0) = 5.0;
        const int change_k = 50;
        const int steps = 150;  // 15 s
        Metrics m = run_scenario(mpc, Ad, Bd, sp_0, sp_1, change_k, steps,
                                 "results_quad_lat_mid.csv");
        print_timing(m);
        check(std::fabs(m.x_final(0, 0) - 5.0) < 0.05, "px settles to 5 m (|error|)",
              std::fabs(m.x_final(0, 0) - 5.0));
        check(m.max_abs_tilt <= TILT_MAX + CONSTR_TOL, "tilt state within limit (max)",
              m.max_abs_tilt);
        check(m.max_abs_tilt_ref <= TILT_MAX + CONSTR_TOL, "tilt ref within limit (max)",
              m.max_abs_tilt_ref);
    }

    // --- Scenario 5: trajectory tracking through the gates ---------------------------------
    std::cout << "=== Scenario 5: trajectory tracking through the gates ===\n";
    {
        MPC mpc(Ad, Bd, C, W_x, W_u, W_du, du_max, du_min, u_max, u_min,
                x_max, x_min, x0, HORIZON, QP);
        const PathPoints points = {{{0.0, 0.0, 0.0},
                                    {3.0, 2.0, 0.0},
                                    {6.0, -2.0, 0.0},
                                    {9.0, 2.0, 0.0},
                                    {12.0, -2.0, 0.0},
                                    {15.0, 0.0, 0.0}}};
        const double speed = 2.0;  // [m/s]
        const Trajectory traj = make_trajectory(points, speed);
        const int steps = 200;     // 20 s
        TrajectoryMetrics tm = run_trajectory(mpc, Ad, Bd, traj, steps,
                                              "results_quad_trajectory.csv");
        const Metrics& m = tm.base;
        print_timing(m);
        std::cout << "  trajectory duration: " << traj.t_arrive[N_PTS - 1U] << " s\n";
        check(tm.max_track_err < 1.0, "tracking error (max ||p - p_ref||)",
              tm.max_track_err);
        for (std::size_t j = 1; j + 1U < N_PTS; ++j) {
            check(tm.min_gate_dist[j] < 0.5, "passes gate (closest approach)",
                  tm.min_gate_dist[j]);
        }
        const double final_err = std::hypot(m.x_final(0, 0) - 15.0,
                                            m.x_final(1, 0), m.x_final(2, 0));
        check(final_err < 0.05, "settles at final point (|error|)", final_err);
        check(m.max_abs_tilt <= TILT_MAX + CONSTR_TOL, "tilt state within limit (max)",
              m.max_abs_tilt);
        check(m.max_abs_tilt_ref <= TILT_MAX + CONSTR_TOL, "tilt ref within limit (max)",
              m.max_abs_tilt_ref);
    }

    std::cout << "\n" << (n_failed == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << " (" << n_failed << " failed)\n";
    return (n_failed == 0) ? 0 : 1;

    
}