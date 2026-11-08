clear; clc;

%% System
A = [0.9002, -0.095;
     0.095,  0.9952];

B = [0.095;
     0.004833];

C = eye(2);
D = zeros(2,1);
Ts = 1;

plant = ss(A, B, C, D, Ts);

%% MPC object
h = 15;
ctrl = mpc(plant, Ts, h, h);

%% Weights
ctrl.Weights.OutputVariables         = [0, sqrt(40)];
ctrl.Weights.ManipulatedVariables    = sqrt(0.1);
ctrl.Weights.ManipulatedVariablesRate = sqrt(0.1);

%% Constraints
ctrl.MV.Min      = -25;
ctrl.MV.Max      =  25;
ctrl.MV.RateMin  = -10;
ctrl.MV.RateMax  =  10;

ctrl.OV(1).Min = -50;
ctrl.OV(1).Max =  50;

ctrl.OV(2).Min = -50;
ctrl.OV(2).Max =  50;

%% Simulation


sim_steps = 150;

x        = zeros(2,1);
xmpc     = mpcstate(ctrl);

setpoint = [0, 10];

% columns:
% k, x_r, x0, x1, u
log = zeros(sim_steps, 5);
ms_tot = 0.0
for k = 1:sim_steps

    if k == 52
        setpoint = [0, 7.5];
    end
    tic
    u = mpcmove(ctrl, xmpc, x', setpoint);
    ms = toc * 1000
    ms_tot = ms_tot + ms
    x = A * x + B * u;

    log(k,:) = [
        k-1, ...
        setpoint(2), ...
        x(1), ...
        x(2), ...
        u
    ];
end

;

fprintf('Time for 1 MPC iteration: %.3f ms\n', ms_tot/sim_steps);

%% Metrics
mse   = mean((log(:,4) - log(:,2)).^2);
u_tot = sum(abs(log(:,5)));

fprintf('MSE  (x1 vs setpoint): %.4f\n', mse);
fprintf('Total control effort: %.4f\n', u_tot);

%% Save CSV for Python
T = array2table(log, ...
    'VariableNames', {'k', 'x_r', 'x0', 'x1', 'u'});

writetable(T, 'results_MATLAB.csv');

fprintf('Saved MATLAB MPC results to results_MATLAB.csv\n');

%% MATLAB Plot
figure('Position', [100 100 900 700]);

subplot(2,1,1);

plot(log(:,1), log(:,3), '--', ...
    'Color', [1 0.5 0], ...
    'DisplayName', 'x0');

hold on;

plot(log(:,1), log(:,4), ...
    'b', ...
    'LineWidth', 2, ...
    'DisplayName', 'x1');

plot(log(:,1), log(:,2), ...
    'r:', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Setpoint');

ylabel('State Values');
title('MATLAB MPC Toolbox');
legend;
grid on;

subplot(2,1,2);

stairs(log(:,1), log(:,5), ...
    'g', ...
    'DisplayName', 'u');

hold on;

yline(25,  'r:', 'DisplayName', 'Max u');
yline(-25, 'm:', 'DisplayName', 'Min u');

xlabel('Time Step (k)');
ylabel('Control Effort');

legend;
grid on;