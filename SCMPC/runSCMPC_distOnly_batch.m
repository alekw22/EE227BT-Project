% Run scenario MPC on the example from Schildbach et al (2014) but with
% only additive disturbances. Run many short simulations.
clear

%% Problem setup
A = [0.7 -0.2; -0.3 0.9];
B = eye(2);
N = 6;          % prediction horizon
K1 = 9;         % number of scenarios for constraint on x1
K2 = 9;         % number of scenarios for constraint on x2
T = 6;          % number of closed loop time steps
M = 200;        % number of close loop simulations
x0 = [1 1]';    % initial condition
sigma = sqrt(0.1);

% allocate space to store closed loop path
xCL = zeros(2,T+1,M);
xCL(:,1,:) = repmat(x0, 1, 1, M);
uCL = zeros(2,T,M);

%% run SCMPC
for m = 1:M
    for t = 1:T
        display([m t])
        % Sample uncertainty
        omega1 = sigma*randn(2,N,K1);
        omega2 = sigma*randn(2,N,K2);
        % Solve problem
        u = solveFHOCP_distOnly(xCL(:,t,m), omega1, omega2);
        % Advance the system
        uCL(:,t,m) = u(:,1);
        xCL(:,t+1,m) = A*xCL(:,t,m) + B*uCL(:,t,m) + sigma*randn(2,1);
    end
end

%% Plot output
figure
hold on
for m = 1:M
    plot(squeeze(xCL(1,:,m)), squeeze(xCL(2,:,m)))
end
plot([1 1], ylim, 'r--')
hold on
plot(xlim, [1 1], 'r--')
xlabel('x1'); ylabel('x2'); title('Closed loop trajectory')

%% Constraint violations
viol = xCL < 1;
numViols = sum(viol,3);
Vemp = numViols/M    % empirical violation rate

%% Closed loop costs
Cmat = zeros(T, M);
for m = 1:M
    for t = 1:T
        Cmat(t,m) = xCL(:,t,m)'*xCL(:,t,m) + uCL(:,t)'*uCL(:,t);
    end
end
Ctot = sum(Cmat)     % total close loop cost
Cavg = mean(Ctot)    % average closed loop stage cost

%% Save
save('data200CL.mat')
