% Run scenario MPC on the example from Cannon et al (2010) on stochastic
% tube MPC. Run many short simulations.
clear

%% Problem setup
A = [1.6 1.1; -0.7 1.2];
B = [1 1]';
g = [1 0.2]';
h = 1.2;
N = 6;           % prediction horizon
K = 4;           % number of scenarios
T = 5;           % number of closed loop time steps
M = 200;         % number of closed loop simulations to run
x0 = [-5 60]';   % initial condition
sigma = 1/12;

% allocate space to store closed loop path
xCL = zeros(2,T+1,M);
xCL(:,1,:) = repmat(x0, 1, 1, M);
uCL = zeros(T,M);

%% run SCMPC
for m = 1:M
    for t = 1:T
        display([m t])
        % Sample uncertainty
        omega = sigma*randn(2,N,K);

        % Solve problem
        [u, x] = solveFHOCP_tubeEx(xCL(:,t,m), omega);

        % Advance the system
        uCL(t,m) = u(1);
        xCL(:,t+1,m) = A*xCL(:,t,m) + B*uCL(t,m) + sigma*randn(2,1);
    end
end

%% Plot output
figure
subplot(121)
hold on
for m = 1:M
    plot(squeeze(xCL(1,:,m)), squeeze(xCL(2,:,m)))
end
x1 = linspace(-5, 2, 50);
x2 = 6 - 5*x1;
plot(x1, x2, 'k')
xlim([-5 2])
xlabel('x1')
ylabel('x2')
title('Closed loop trajectories')

subplot(122)
hold on
for m = 1:M
    plot(squeeze(xCL(1,:,m)), squeeze(xCL(2,:,m)), '.')
end
x1 = linspace(-5, 2, 50);
x2 = 6 - 5*x1;
plot(x1, x2, 'k')
xlim([-2.5 -1.5])
ylim([13 19])
xlabel('x1')
ylabel('x2')
title('Close up: First step')

%% Analyze results
% violations
numViols = zeros(1,T+1);
for t = 1:(T+1)
    viol = g'*squeeze(xCL(:,t,:)) > h;
    numViols(t) = sum(viol>0);
end
Vemp = numViols/M    % empirical violation rate

% stage costs
Cmat = zeros(T,M);
for t = 1:T
    for m = 1:M
        Cmat(t,m) = xCL(:,t,m)'*xCL(:,t,m) + uCL(t,m)'*uCL(t,m);
    end
end
Cavg = mean(Cmat,2)   % average closed loop stage cost
sum(Cavg)

%% Save
% save('data200_tubeEx.mat')
