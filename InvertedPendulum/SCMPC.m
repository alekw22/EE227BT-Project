% Do scenario MPC on the inverted pendulum system
clear all

%% Setup
load('sysDynDamped.mat')
% Q = 10*eye(4);
Q = [1 0 0 0;
     0 0 0 0;
     0 0 1 0;
     0 0 0 0];
R = 0.05;
sigma = 0.05;
K = 2;
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
end

%%
omega1 = sigma*randn(N,K);
xCL = zeros(4,N+1,K);
xCL(:,1,:) = repmat(x0, 1, 1, K);

for n = 1:N
    xCL(:,n+1,1) = A*xCL(:,n) + B*omega1(n,1)
    xCL(:,n+1,2) = A*xCL(:,n) + B*omega1(n,2)
end

%%
figure;
plot(xCL(1,:,1))
hold on
plot(xCL(1,:,2))

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

%%
figure
plot((0:T)*0.05, xCL(1,:))
hold on
for t = 1:5:T
    plot(((0:N)+t-1)*0.05, xOL(1,:,t), 'k--')
end

%%
figure
plot((0:T)*0.02, xCL(3,:), 'bo')
hold on
for t = 1:5:T
    plot(((0:N)+t-1)*0.05, xOL(3,:,t), 'k--')
end

