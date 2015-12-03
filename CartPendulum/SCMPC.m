% Do scenario MPC on the cart pendulum system
clear all

%% Setup
load('sysDyn.mat')
% Q = 10*eye(4);
Q = [1 0 0 0;
     0 0 0 0;
     0 0 1 0;
     0 0 0 0];
R = 0.05;
sigma = 0.05;
K = 4;
N = 20;
T = 150;
x0 = [1 0 0 0]'; % start the cart away from the origin

% allocate space to store closed loop path
xCL = NaN(4,T+1);
xCL(:,1) = x0;
uCL = NaN(1,T);

%% Solve
for t = 1:T
    t
    omega = sigma*randn(N,K);
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
            x(3,2:end,:) <= 0.1;
            -x(3,2:end,:) <= 0.1;
            for k = 1:K
                x(:,1,k) == xCL(:,t);
                for n = 1:N
                    x(:,n+1,k) == A*x(:,n,k) + B*(u(n) + omega(n,k));
                end
            end
    cvx_end
    
    % store closed loop results
    uCL(t) = u(1);
    xCL(:,t+1) = A*xCL(:,t) + B*(u(1) + sigma*randn);
    display(xCL(:,t+1)')
end

%% Plot
figure
subplot(211)
plot((0:T)*0.05, xCL(1,:))
ylabel('x')

subplot(212)
plot((0:T)*0.05, xCL(3,:))
ylabel('\phi')
xlabel('Time')

figure
plot((0:(T-1))*0.05, uCL)
ylabel('Input')
xlabel('Time')

%% Violations
numViol = sum(abs(xCL(3,:))>0.1)
