% Do scenario MPC on the cart pendulum system
clear all

%% Setup
load('sysDyn.mat')
Q = [1 0 0 0;
     0 0 0 0;
     0 0 1 0;
     0 0 0 0];
R = 0.05;
wmax = 0.15;
K = 4;
N = 6;
T = 150;
M = 100;
x0 = [1 0 0 0]'; % start the cart away from the origin

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
                x(3,2:end,:) <= 0.1;
                -x(3,2:end,:) <= 0.1;
                for k = 1:K
                    x(:,1,k) == xCL(:,t,m);
                    for n = 1:N
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
plot((0:T)*0.05, xCL(1,:,m))
ylabel('x')

subplot(212)
plot((0:T)*0.05, xCL(3,:,m))
ylabel('\phi')
xlabel('Time')

figure
plot((0:(T-1))*0.05, uCL(:,m))
ylabel('Input')
xlabel('Time')

%% Violations
numViol = sum(abs(xCL(3,:))>0.1)
%%
save('resultsN6.mat')
clear xCL uCL

% allocate space to store closed loop path
N = 15;
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
                x(3,2:end,:) <= 0.1;
                -x(3,2:end,:) <= 0.1;
                for k = 1:K
                    x(:,1,k) == xCL(:,t,m);
                    for n = 1:N
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

%%
save('resultsN15.mat')
