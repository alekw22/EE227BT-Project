% Do scenario MPC on the cart pendulum system
clear all

%% Setup
load('sysDyn.mat')     % load A, B matrices
Q = [1 0 0 0;          % weight on the state in the cost function
     0 0 0 0;
     0 0 1 0;
     0 0 0 0];
R = 0.05;              % weight on the input
wmax = 0.15;           % maximum disturbance magnitude
K = 4;                 % sample complexity
N = 6;                 % prediction horizon
T = 150;               % closed-loop time steps
M = 50;                % number of systems to simulate
x0 = [1 0 0 0]';       % start the cart away from the origin

% allocate space to store closed loop path
xCL = NaN(4,T+1,M);
xCL(:,1,:) = repmat(x0, 1, 1, M);
uCL = NaN(T,M);

%% Solve
for m = 1:M
    for t = 1:T
        display([N m t])
        omega = 2*wmax*(rand(N,K)-0.5);
        cvx_begin
            variables x(4,N+1,K) u(1,N)
            objfun = 0;
            for k = 1:K
                for n = 1:N
                    objfun = objfun + x(:,n,k)'*Q*x(:,n,k) + u(n)*R*u(n);
                end
            end
            minimize(objfun)
            subject to
                x(3,2:end,:) <= 0.1;     % state constraint 1
                -x(3,2:end,:) <= 0.1;    % state constraint 2
                for k = 1:K
                    x(:,1,k) == xCL(:,t,m);    % initial condition
                    for n = 1:N
                        % system dynamics
                        x(:,n+1,k) == A*x(:,n,k) + B*(u(n) + omega(n,k));
                    end
                end
        cvx_end

        % store closed loop results
        uCL(t,m) = u(1);
        xCL(:,t+1,m) = A*xCL(:,t,m) + B*(u(1) + 2*wmax*(rand()-0.5));
        display(xCL(:,t+1,m)')
    end
end

%% Plot
figure
subplot(211)
plot((0:T)*0.05, squeeze(xCL(1,:,:)))
title('Closed-loop trajectory')
ylabel('x [m]')

subplot(212)
plot((0:T)*0.05, squeeze(xCL(3,:,:)))
hold on
plot(xlim, [0.1 0.1; -0.1 -0.1]', 'k')
ylabel('\theta [rad]')
xlabel('Time [s]')

figure
plot((0:(T-1))*0.05, uCL)
title('Input')
ylabel('Force [N]')
xlabel('Time [s]')

%% Violations
viols = abs(squeeze(xCL(3,:,1:50))) > 0.1;
numViol = sum(viols,2);
Vemp = numViol/M;

figure
plot((0:T)*0.05, Vemp)
ylabel('Percent')
xlabel('Time [s]')
title('Empirical Violation Level')

%% Cost
Cmat = zeros(T,M);
for m = 1:M
    for t = 1:T
        Cmat(t,m) = xCL(:,t,m)'*Q*xCL(:,t,m) + R*uCL(t,m)^2;
    end
end

Cavg = mean(Cmat,2);
Ctot = sum(Cavg)


%%
% save('resultsN6.mat')
