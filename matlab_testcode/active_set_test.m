%% compare_active_set.m
%
% Solves the same 4 QP test problems as tests/test_active_set.cpp using
% MATLAB's quadprog, so the results can be diffed against the C++ output.
%
% quadprog convention:  min 0.5*x'*H*x + f'*x   s.t.  A*x <= b
% This matches the Active_Set class exactly (Q -> H, c -> f).

clear; clc;

opts = optimoptions('quadprog', 'Display', 'off');

%% Test 1: Interior optimum, no active constraints
fprintf('=== Test 1: Interior optimum, no active constraints ===\n');
Q1 = [2 0; 0 2];
c1 = [-2; -2];
A1 = [1 1; -1 0; 0 -1];
b1 = [10; 0; 0];
[x1, fval1] = quadprog(Q1, c1, A1, b1, [], [], [], [], [], opts);
disp('x ='); disp(x1');
fprintf('objective = %.6f\n', fval1);
fprintf('expected:  x = [1, 1], objective = -2\n\n');

%% Test 2: Single active constraint (projection)
fprintf('=== Test 2: Single active constraint (projection) ===\n');
Q2 = [1 0; 0 1];
c2 = [-1; -1];
A2 = [1 1; -1 0; 0 -1];
b2 = [1; 0; 0];
[x2, fval2] = quadprog(Q2, c2, A2, b2, [], [], [], [], [], opts);
disp('x ='); disp(x2');
fprintf('objective = %.6f\n', fval2);
fprintf('expected:  x = [0.5, 0.5], objective = -0.75\n\n');

%% Test 3: Nocedal & Wright classic QP
fprintf('=== Test 3: Nocedal & Wright classic QP ===\n');
Q3 = [2 0; 0 2];
c3 = [-2; -5];
A3 = [-1 2; 1 2; 1 -2; -1 0; 0 -1];
b3 = [2; 6; 2; 0; 0];
[x3, fval3] = quadprog(Q3, c3, A3, b3, [], [], [], [], [], opts);
disp('x ='); disp(x3');
fprintf('internal objective = %.6f\n', fval3);
fprintf('f(x) = (x1-1)^2+(x2-2.5)^2 = %.6f\n', (x3(1)-1)^2 + (x3(2)-2.5)^2);
fprintf('expected:  x = [1.4, 1.7], f(x) = 0.8\n\n');

%% Test 4: 3D QP with coupled variables
fprintf('=== Test 4: 3D QP with coupled variables ===\n');
Q4 = [4 1 0; 1 2 0; 0 0 2];
c4 = [-6; -3; -4];
A4 = [1 1 1; -1 0 0; 0 -1 0; 0 0 -1];
b4 = [3; 0; 0; 0];
[x4, fval4] = quadprog(Q4, c4, A4, b4, [], [], [], [], [], opts);
disp('x ='); disp(x4');
fprintf('objective = %.6f\n', fval4);

fprintf('Active constraints (residual ~0 means active):\n');
disp(A4 * x4 - b4);
