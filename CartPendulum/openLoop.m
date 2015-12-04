% Open loop simulation to check that the dynamics make sense
clear
load('sysDyn.mat')

%% Simulate autonomous system
x0 = [0 0 0.2 0]';
T = 200;
x = [x0, zeros(4,T)];

for t = 1:T
    x(:,t+1) = A*x(:,t);
end

%% Plot
figure
subplot(211)
plot((0:T)*0.05, x(1,:))
ylabel('x')

subplot(212)
plot((0:T)*0.05, x(3,:))
ylabel('\theta')
xlabel('Time')

%% Simulate system subject to disturbances
K = 50;   % number of systems
T = 20;   % number of time steps
x0 = [0 0 0 0]';
x = zeros(4,T+1,K);
x(:,1,:) = repmat(x0, 1, K);
sigma = 0.25;

for t = 1:T
    for k = 1:K
        x(:,t+1,k) = A*x(:,t,k) + Bw*sigma*randn();
    end
end

%% Plot
figure
subplot(211)
plot((0:T)*0.05, squeeze(x(1,:,:)))
ylabel('x')

subplot(212)
plot((0:T)*0.05, squeeze(x(3,:,:)))
ylabel('\theta')
xlabel('Time')

