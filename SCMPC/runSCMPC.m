% Run scenario MPC on the example from Schildbach et al (2014) that has
% both additive and multiplicative disturbances
clear

%% Problem setup
Anom = [0.7 -0.2; -0.3 0.9];   % nominal state transition matrix
Aunc = [0 -0.1; -0.2 0];       % coefficients of random variable theta
B = eye(2);
N = 5;         % prediction horizon
K = 19;        % number of scenarios
T = 2500;      % number of closed loop time steps
x0 = [1 1]';   % initial condition
sigma = 0.1;

% allocate space to store closed loop path
xCL = zeros(2,T+1);
xCL(:,1) = x0;
uCL = zeros(2,T);

%% run SCMPC
for t = 1:T
    t
    datetime()
    % Sample uncertainty
    theta = rand(N,K);
    omega = sigma*randn(2,N,K);

    % Solve problem
    [u, x] = solveFHOCP(xCL(:,t), theta, omega);
    
    % Advance the system
    uCL(:,t) = u(:,1);
    xCL(:,t+1) = (Anom + Aunc*rand())*xCL(:,t) + B*uCL(:,t) + sigma*randn(2,1);
end

%% Plot output
figure
plot(xCL(1,:), xCL(2,:), '.')
hold on
plot([1 1], ylim, 'r--')
hold on
plot(xlim, [1 1], 'r--')
xlabel('x1'); ylabel('x2'); title('Closed loop trajectory')

%% Analyze results
% violations
viol = xCL < 1;
numViols = sum(sum(viol)>0);
Vemp = numViols/T    % empirical violation rate

% stage costs
Cvect = zeros(1,T);
for t = 1:T
    Cvect(t) = xCL(:,t)'*xCL(:,t) + uCL(:,t)'*uCL(:,t);
end
Cavg = mean(Cvect)   % average closed loop stage cost
Cstd = std(Cvect)    % std. deviation of closed loop stage cost

%% Save
% save('data2500_2.mat')
