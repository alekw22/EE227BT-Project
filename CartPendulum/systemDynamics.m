% Discrete time dynamics of a pendulum attached to a cart
clear

% System Parameters
M = 0.5;       % mass of cart [kg]
m = 0.2;       % mass of the pendulum [kg]
b = 0.1;       % cart coefficient of friction [N/m/s]
c = 0.01;      % pendulum coefficient of friction [N.m.s]
I = 0.006;     % pendulum moment of inertia [kg.m/s^2]
g = 9.8;       % acceleration of gravity [m/s^2]
l = 0.3;       % half length of pendulum [m]
T = 0.05;      % discretization time [s]

% Continuous time dynamics
d = I*(M+m) + M*m*l^2;
c1 = -(I+m*l^2)*b/d;
c2 = ((m^2)*g*l^2)/d;
c3 = c*m*l/d;
c4 = m*l*b/d;
c5 = -m*g*l*(M+m)/d;
c6 = -c*(M+m)/d;
c7 = (I+m*l^2)/d;
c8 = -m*l/d;

Ac = [0 1 0 0;
      0 c1 c2 c3;
      0 0 0 1;
      0 c4 c5 c6];
Bc = [0 c7 0 c8]';
Cc = eye(4);
Dc = zeros(4,1);

ctSys = ss(Ac,Bc,Cc,Dc);

% Discrete time dynamics, calculated using zero order hold
sys = c2d(ctSys, T, 'zoh');
A = sys.A;
B = sys.B;
Bw = B;

%%
% save('sysDyn.mat', 'A', 'B', 'Bw')
