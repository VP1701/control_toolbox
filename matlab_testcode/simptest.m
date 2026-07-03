% ==================================================================
%  MATLAB REFERENCE SOLVER – 6 test problems for your Simplex code
%  Run this script and compare with your C++ implementation
% ==================================================================
clear; clc; close all;

problems = {
    % 1. Textbook Example
    struct('name','Textbook example (max 2x1+3x2, x1+x2<=4, 2x1+3x2<=9)', ...
           'A',[2 3; -1 1], 'b',[6;1], 'c',[-1;-3], 'expected',[], 'z',[]), ...
    
    % 2. Infeasible
    struct('name','Infeasible problem', ...
           'A',[1 1; 1 -1], 'b',[1; -2], 'c',[-1;-1], 'expected',[], 'z',[]), ...
    
    % 3. Unbounded
    struct('name','Unbounded problem', ...
           'A',[-1 -1], 'b',[-1], 'c',[-1;-2], 'expected',[], 'z',-Inf), ...
    
    % 4. Degeneracy
    struct('name','Degeneracy test (multiple optimal solutions)', ...
           'A',[1 4; 1 2], 'b',[8;4], 'c',[-3;-9], 'expected',[], 'z',-3), ...
    
    % 5. Large-ish Random Problem (Fixed A matrix and b vector)
    struct('name','Large-ish random feasible problem (10 vars, 8 constraints)', ...
           'A',[ ...
            1.484, 0.937, 2.744, 1.095, 0.632, 2.173, 1.918, 0.049, 2.501, 1.723; ...
            0.294, 2.181, 1.639, 0.917, 2.884, 0.425, 1.273, 2.660, 0.167, 1.984; ...
            2.463, 0.751, 1.852, 2.094, 0.388, 1.629, 2.987, 0.883, 1.526, 0.279; ...
            1.172, 2.695, 0.604, 1.483, 2.350, 0.976, 0.251, 1.840, 2.712, 0.668; ...
            0.895, 1.329, 2.108, 0.447, 1.763, 2.584, 0.732, 1.195, 0.339, 2.876; ...
            2.017, 0.583, 1.964, 2.331, 0.108, 1.472, 2.695, 0.624, 1.883, 0.950; ...
            1.706, 2.249, 0.417, 1.628, 2.973, 0.862, 1.354, 2.108, 0.775, 1.491; ...
            0.628, 1.883, 2.540, 0.296, 1.017, 2.762, 0.194, 1.650, 2.429, 0.583 ...
           ], ...
           'b',[18.743; 19.294; 19.826; 18.472; 19.105; 19.683; 20.351; 19.628], ... % Converted to column vector (;)
           'c',[-2.351, -3.827, -1.194, -3.506, -2.683, -1.927, -3.148, -2.074, -1.583, -3.992], ...
           'expected',[], 'z',[]), ...
           
    % 6. Refinery Blending Problem (Feasible Version)
    struct('name','Refinery blending optimization problem (Feasible)', ...
       'A', [ ...
        % Feedstock (<=)
         1,  1,  1,  0,  0,  0,  0,  0,  0; ... % Low-tier <= 20000
         0,  0,  0,  1,  1,  1,  0,  0,  0; ... % Mid-tier <= 20000
         0,  0,  0,  0,  0,  0,  1,  1,  1; ... % High-tier <= 10000
        % Market Demand (>= 20000, 10000, 5000 -> <= -20000, -10000, -5000)
        -1,  0,  0, -1,  0,  0, -1,  0,  0; ... % Regular
         0, -1,  0,  0, -1,  0,  0, -1,  0; ... % Mid-grade
         0,  0, -1,  0,  0, -1,  0,  0, -1; ... % Premium
        % Octane (>= 0 -> <= 0)
         3,  0,  0, -2,  0,  0, -6,  0,  0; ... % Regular
         0,  5,  0,  0,  0,  0,  0, -4,  0; ... % Mid-grade
         0,  0,  9,  0,  4,  0,  0,  0,  0; ... % Premium
        % Sulfur (<= 0)
         0.2, 0, 0, -0.6, 0, 0, -0.9, 0, 0; ... % Regular
         0, 0.5, 0, 0, -0.3, 0, 0, -0.6, 0; ... % Mid-grade
         0, 0, 0.8, 0, 0, 0, 0, 0, -0.3     ... % Premium
       ], ...
       'b', [20000; 20000; 10000; -20000; -10000; -5000; 0; 0; 0; 0; 0; 0], ...
       'c', -[125-58, 132-58, 140-58, 125-65, 132-65, 140-65, 125-74, 132-74, 140-74], ...
       'expected',[], 'z',[])
};

fprintf('=== MATLAB REFERENCE SOLUTIONS ===\n\n');
for i = 1:length(problems)
    prob = problems{i};
    
    fprintf('Problem %d: %s\n', i, prob.name);
    
    c_vec = prob.c(:);
    
    options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
    [xopt, fopt, exitflag, output] = linprog(c_vec, prob.A, prob.b, [], [], ...
                                             zeros(length(c_vec),1), [], options);
    
    if exitflag == 1
        fprintf('MATLAB says: OPTIMAL\n');
        fprintf('x* = ['); fprintf(' %.4f', xopt); fprintf(' ]\n');
        fprintf('objective = %.10g\n', fopt);
    elseif exitflag == -2
        fprintf('MATLAB says: INFEASIBLE\n');
    elseif exitflag == -3
        fprintf('MATLAB says: UNBOUNDED\n');
    else
        fprintf('MATLAB exitflag = %d (%s)\n', exitflag, output.message);
    end
    
    
end