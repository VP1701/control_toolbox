# Control Toolbox

This is my hobby repository for doing control theory related C++ programs.



## Currently working on
- extending matrix class
- LTV-MPC
- adding CMake

## future plans
- State feedback controller
- State predictor
- Kitagawa algorithm
- Kalman predictor
- Kalman filter

## MPC solver for Linear time invariant MIMO systems
compile code:
```g++ -O2 -Iinclude tests/mpc_test.cpp src/simplex.cpp src/matrix.cpp src/active_set.cpp src/mpc.cpp -o test_mpc```
run the code:
```./test_mpc```

## Active set algorithm
compile code:
```g++ -O2 -Iinclude tests/test_active_set.cpp src/simplex.cpp src/matrix.cpp src/active_set.cpp -o test_active_set```
run the code:
```./test_active_set```


## Simplex algorithm


For testing the solver you need to clone the repo, have C++ and g++ installed. After that you can run the test with these commands:

compile code:
```g++ -O2 -Iinclude tests/test_simplex.cpp src/simplex.cpp src/matrix.cpp -o test_simplex```
run the code:
```./test_simplex```


