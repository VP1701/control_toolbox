# Control Toolbox (C++)

This is my hobby repository for doing control theory related C++ programs.

## Key Features
- Model Predictive Control (MPC): Formulated for Linear TIme-invariant (LTI) MIMO systems. Supports Linear and Quadratic cost functions.
- Custom optimization algorithms
    - Simplex algorithm:
    - Active Set algortihm: 
- PID control: Both algorithm support antiwindup
    - Positional PID-control algorithm
    - Velocity form PID-control algorithm

## Validation and results
The C++ implementations in this library are verifies against MATLAB. Simplex- and Active Set -algortihms were tested against matlab
and gave same results. The MPC-algortihm was also tested against matlab and gave almost identical results. 

MPC-algortihm was tested on a 2 state integrator system and as seen below the results are almost identical to matlabs with the QP-cost
![MPC Comparison](mpc_comparison.png)

## Building tests

This project uses CMake. From the repo root:

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

This builds `simplex_test`, `active_set_test`, and `mpc_test` into the `build/` folder. 

## Running tests
### MPC solver for Linear time invariant MIMO systems
run the code:
```
./mpc_test
```

### Active set algorithm
run the code:
```
./active_set_test
```

### Simplex algorithm
run the code:
```
./simplex_test
```


## Roadmap and future plans
I'm Currently working on implementing a system/plant class to streamline simulations and system definitions. Improvements and test for the PID class are also being worked on.

### Next features
- State feedback controller
- State predictor
- Kitagawa algorithm
- Kalman predictor
- Kalman filter