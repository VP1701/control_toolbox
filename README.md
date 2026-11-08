# Control Toolbox (C++)

This is my personal passion project C++ library for modern control and optimization algorithms. The goal of this repository is to learn
control theory and optimization and the implementation of these in C++. I am building the algorithms from scratch myself to learn how they work.
The rule for this project is to never use any other external libraries than the C++ standard library. The goal is that the developed algorithms are
robust and efficient enough to be used in real-world applications.

## Key Features
- Model Predictive Control (MPC): Formulated for Linear Time-invariant (LTI) MIMO systems. Supports Linear and Quadratic cost functions.
- Custom optimization algorithms
    - Simplex algorithm
    - Active Set algorithm
- PID control: Both algorithms support antiwindup
    - Positional PID-control algorithm
    - Velocity form PID-control algorithm
- Kalman filter

## Validation and results
Simplex- and Active Set -algorithms were tested against MATLAB
and gave same results. The MPC algortihm was also tested against MATLAB and gave almost identical results. 

MPC algortihm was tested on a stable 2-state system and as seen below the results are almost identical to MATLAB's with the QP-cost
![MPC Comparison](documentation/images/mpc_comparison.png)


After MPC validation against MATLAB the MPC has been tested on a 8 state drone model which was linearized around hover. 
states:
x = [px, py, pz, vx, vy, vz, phi, theta]^T
controls:
u = [phi_r, theta_r, T]^T
where phi is the roll and theta is the pitch and phi_r and theta_r are the setpoints given to the attitude controller. The attitude controller is part of the model and it is modeled as a second order system with dead time approximated as a first order system. Source for state space model is "Model Predictive Control for Trajectory Tracking of Unmanned Aerial Vehicles Using Robot Operating System" by Mina Kamel, Thomas Stastny, Kostas Alexis and Roland Siegwart (2017). The source had an error in printed out matrix where attitude columns were swapped.

Parameters used are illustrative and not identified from a real drone and the state space system was discretized with MATLAB's c2d using ZOH.

The MPC had horizon of 20 steps (2 s), tilt was limited to +-20 degrees, tilt rate to 2 rad/s and thrust between 0.3g and 1.8 g. The trajectory was designed to be flown at 2 m/s. Results on the simulation can be seen below.

![MPC on drone trajectory](documentation/images/linear_drone_trajectory_track.png)

The tests were done on the exact same plant which the MPC used and no disturbances were added. In reality there would be model mismatch and noise. The trajectory starts at the full 2 m/s and drone starts at 0 m/s which causes some deviation from trajectory at the first segment before first gate. The solve time stays under 10 ms with worst case solve being 6.1 ms on my laptop with Ryzen 5 4500U. In the future this will also be tested on the Jetson Orin Nano with the same goal of worst case solve time being under 10 ms.

## Building tests (on Linux)

This project uses CMake. From the repo root:

```bash
mkdir build && cd build
cmake ..
cmake --build .
cd ..
```

This builds the different test cases in the project into the `build/` folder. 


Create a Python virtual environment for the plotters
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running tests

All commands are run from the repo root. Activate the Python environment
first (`source .venv/bin/activate`) to use the plotters. The tests write
their CSV files to the folder they are run from.

### MPC on a 2-state system, compared against MATLAB
Runs the MPC with both the linear (LP) and quadratic (QP) cost. The plotter
compares both against MATLAB's MPC Toolbox reference
(`matlab_testcode/results_MATLAB.csv`) and prints the tracking error and
control effort of each.
```bash
./build/mpc_test
python3 plot_mpc_test.py
```

### Active set QP solver
Solves a test QP and prints the result next to MATLAB's `quadprog` solution.
```bash
./build/active_set_test
```

### Simplex LP solver
Solves a test LP and prints the result next to MATLAB's solution.
```bash
./build/simplex_test
```

### Kalman filter
Estimates position and velocity from noisy position measurements. The plot
is saved to `kalman_results.png`.
```bash
./build/kalman_demo
python3 plot_kalman.py
```

### PID controllers
Runs the positional and velocity form PID controllers.
```bash
./build/pid_test
python3 plot_pid_test.py
```

### Hover-linearized drone (MPC)
Runs five scenarios on the drone model: hover hold, altitude step, lateral step,
mid-run lateral step and trajectory tracking through gates. Prints PASS/FAIL
for each check and the solve times.
```bash
./build/linear_drone
python3 plot_quad_test.py              # scenarios 1-4, add --save for PNGs
python3 plot_quad_trajectory.py        # trajectory tracking, add --save for a PNG
```

## Roadmap and future plans
### I'm currently working on:
- Implementing a system/plant class to streamline simulations and system definitions. 
- Kitagawa algorithm
- PID tuning in pid_test.cpp for fair comparison with MPC

### Next steps:
- Optimizing the active set algorithm. KKT factorization could be updated incrementally instead of doing full factorization every iteration.
- Status handling for the MPC and matrix classes

### Drone MPC:
- Timing on Jetson Orin Nano
- Adding disturbances to simulation and implementing state estimation with Kalman filter

### Long-term goals:
- Making separate STM32 compatible code and modifying matrix class and other algorithms to work on STM32H753.
- Testing on a NUCLEO-H753ZI board

