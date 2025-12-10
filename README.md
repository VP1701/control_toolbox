# Control Toolbox

This is my hobby repository for doing control theory related C++ scripts and/or programs.


## MATH.720 Modeling report
my simplex algorithm implementation is in this repository. The simplex algorithm minimizes the given problem and calculates the optimal objective value. The implementation also handles degenerate cases/cycling with Bland's rule. The solver isn't perfect the GEQ constraint type doesn't work properly so you still have to manually flip the inequalities to LEQ which doesn't prevent the use of the algorithm but makes it minutely more inconvenient to use.

Inside main.cpp there are test cases for the solver. Last one of them is the soltuion to the gasoline refining example from my report. Because the simplex solver I made is made for minimizing an objective function the solver returns negative of the maximum profits. but the value given by the solver * -1 is the correct answer. The solvers solution were also checked with matlabs linprog and it gave the same answer so I can confidently say that my solver works.

For testing the solver you need to clone the repo, have C++ and g++ installed. After that you can run the test with these commands:

you can run it with these commands
```g++ -Iinclude src/main.cpp src/simplex.cpp src/matrix.cpp -o test```
```./test.exe```


## done
- Discrete PID-controller class that includes different version depending on the chosen discretization method.
    includes antiwindup if chosen

## Currently working on
- extending matrix class
- start working on vector class

## future plans
-State feedback controller
-State predictor
-Kalman predictor
-Kalman filter