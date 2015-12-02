clear

%% Problem setup
A = [0.7 -0.2; -0.3 0.9];
B = eye(2);
N = 6;           % prediction horizon
K1 = 9;          % number of scenarios for constraint on x1
K2 = 9;          % number of scenarios for constraint on x2
T = 3000;        % number of closed loop time steps
x0 = [1 1]';     % initial condition
sigma = sqrt(0.1);

% allocate space to store closed loop path
xCL = zeros(2,T+1);
xCL(:,1) = x0;
uCL = zeros(2,T);

%% run SCMPC
for t = 1:T
    t

    % Sample uncertainty
    omega1 = sigma*randn(2,N,K1);
    omega2 = sigma*randn(2,N,K2);

    % Solve problem
    u = solveFHOCP_distOnly(xCL(:,t), omega1, omega2);
    
    % Advance the system
    uCL(:,t) = u(:,1);
    xCL(:,t+1) = A*xCL(:,t) + B*uCL(:,t) + sigma*randn(2,1);
end

%% Plot output
figure; hold on
plot(xCL(1,:), xCL(2,:), '.')
xlim([0 4])
ylim([0 4])
plot([1 1], ylim, 'k')
hold on
plot(xlim, [1 1], 'k')
% xlabel('x1'); ylabel('x2'); title('Closed loop trajectory')

%% Constraint violations
viol = xCL < 1;
numViols = sum(viol,2);
Vemp = numViols/T    % empirical violation rate

%% Closed loop cost
Cvect = zeros(1,T);
for t = 1:T
    Cvect(t) = xCL(:,t)'*xCL(:,t) + uCL(:,t)'*uCL(:,t);
end
Ctot = sum(Cvect)    % total close loop cost
Cavg = mean(Cvect)   % average closed loop stage cost
Cstd = std(Cvect)    % std. deviation of closed loop stage cost

%% Save
% save('data2500_2.mat')
