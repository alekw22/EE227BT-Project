% Solve a single step of the SCMPC example problem from Schildbach et al.
% (2014). Plot all of the open loop trajectories for each sampled scenario.
clear

%% Sample uncertainty
N = 5;
K = 19;
theta = rand(N,K);
omega = 0.1*randn(2,N,K);
x0 = [2; 2];

%% Solve problem
[u, x] = solveFHOCP(x0, theta, omega);

%% Plot output
figure
for k = 1:K
    plot(squeeze(x(1,:,k)), squeeze(x(2,:,k)))
    hold on
end
xlabel('x1'); ylabel('x2'); title('Open loop trajectories')

figure
plot(0:(N-1), u')
xlabel('Time'); ylabel('Input')
legend('u1', 'u2')