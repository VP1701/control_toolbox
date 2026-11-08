% linear_drone_discretization.m

% source for the drone state space matrices is from ....
% The parameter values are illustrative
% The attitude control is modeled as a second order system. The attitude
% dynamics are assumed to be identical in pitch and roll.

clear all

% drone parameters
g = 9.81; % [m/s^2] gravity
a = 0.01; % [1/s] mass-normalized linear drag. 
% assumed identical on all axes
A_x = a;
A_y = a;
A_z = a;

K_att = 1.0; % [-] DC gain of attitude loop
zeta = 0.7; % [-] attitude loop damping
wn = 12.0; % [rad/s] attitude loop natural frequency
Td = 0.01 % [s] attitude dead time
tau = 2 * zeta/wn + Td; % [s] first order equivalent time constant
Ts = 0.1; % [s] MPC sample time

K_theta = K_att;
K_phi = K_att;
tau_theta = tau;
tau_phi = tau;

% continuous model
nx = 8; % amount of states
nu = 3; % amount of controls
A = zeros(nx);
A(1,4) = 1.0;
A(2,5) = 1.0;
A(3,6) = 1.0;
% In the source the printed matrix has these two swapped
A(4,8) = g;
A(5,7) = -g;

A(4,4) = -A_x;
A(5,5) = -A_y;
A(6,6) = -A_z;
A(7,7) = -1/tau_phi;
A(8,8) = -1/tau_theta;

B = zeros(nx, nu);
B(6,3) = 1.0;
B(7,1) = K_phi/tau_phi;
B(8,2) = K_theta/tau_theta;

% all of the states are assumed to be observable
C = eye(nx);

D = zeros(nx, nu)
% state space system is discretized with zoh assumption
sys = ss(A,B,C,D);
sys_d = c2d(sys,Ts, 'zoh') % these are copied to the cpp test script