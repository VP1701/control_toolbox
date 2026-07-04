# Control Toolbox

This is my hobby repository for doing control theory related C++ programs.

## Currently working on
- system/plant class for making simulation easier

## future plans
- State feedback controller
- State predictor
- Kitagawa algorithm
- Kalman predictor
- Kalman filter

## Building tests

This project uses CMake. From the repo root:

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

This builds `simplex_test`, `active_set_test`, and `mpc_test` into the `build/` folder. 

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