% Discrete time system dynamics

% System Parameters
M = 0.5;       % mass of cart [kg]
m = 0.2;       % mass of the pendulum [kg]
b = 0.1;       % cart coefficient of friction [N/m/s]
I = 0.006;     % pendulum moment of inertia [kg.m/s^2]
g = 9.8;       % acceleration of gravity [m/s^2]
l = 0.3;       % length of pendulum
T = 0.02;      % discretization time

% Discrete time dynamics
d = I*(M+m) + M*m*l^2;
c1 = -(I+ml^2)*b/d;
c2 = ((m^2)*g*l^2)/d;
c3 = -m*l*b/d;
c4 = m*g*l*(M+m)/d;
c5 = I+m*l^2/d;
c6 = m*l/d;

A = [1 T 0 0;
    0 1+c1*T c2*T 0;
    0 0 1 T;
    0 c3*T c4*T 1];

B = [0 c5*T 0 c6*T]';
Bw = B;
