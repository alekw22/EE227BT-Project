% Discrete time system dynamics

% System Parameters
M = 0.5;       % mass of cart [kg]
m = 0.2;       % mass of the pendulum [kg]
b = 0.1;       % cart coefficient of friction [N/m/s]
I = 0.006;     % pendulum moment of inertia [kg.m/s^2]
g = 9.8;       % acceleration of gravity [m/s^2]
l = 0.3;       % length of pendulum
T = 0.02;      % discretization time

% Continuous time dynamics
d = I*(M+m) + M*m*l^2;
c1 = -(I+m*l^2)*b/d;
c2 = ((m^2)*g*l^2)/d;
c3 = -m*l*b/d;
c4 = m*g*l*(M+m)/d;
c5 = (I+m*l^2)/d;
c6 = m*l/d;

Ac = [0 1 0 0;
      0 c1 c2 0;
      0 0 0 1;
      0 c3 c4 0];
Bc = [0 c5 0 c6]';
Cc = eye(4);
Dc = zeros(4,1);

ctSys = ss(Ac,Bc,Cc,Dc);

% Discrete time dynamics
sys = c2d(ctSys, T, 'zoh');
A = sys.A;
B = sys.B;
Bw = B;

%%
% save('sysDyn.mat', 'A', 'B', 'Bw')
