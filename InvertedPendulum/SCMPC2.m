% Do a deterministic MPC on the inverted pendulum system
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
N = 20;
T = 150;
x0 = [1 0 0 0]'; % start the cart away from the origin

% allocate space to store closed loop path
xOL = NaN(4,N+1,T);
xCL = NaN(4,T+1);
xCL(:,1) = x0;
uCL = NaN(1,T);

%% Solve
for t = 1:T
    t
    cvx_begin
        variables x(4,N+1) u(1,N)
        objfun = 0;
        for n = 1:N
            objfun = objfun + x(:,n)'*Q*x(:,n) + u(n)*R*u(n);
        end
        minimize(objfun)
        subject to
        x(:,1) == xCL(:,t);
%         x(3,:) <= 0.10;
%         -x(3,:) <= 0.10;
        for n = 1:N
            x(:,n+1) == A*x(:,n) + B*u(n);
        end
    cvx_end
    
    % store closed loop results
    uCL(t) = u(1);
    xCL(:,t+1) = A*xCL(:,t) + B*(u(1) + sigma*randn);
    xOL(:,:,t) = x;
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

