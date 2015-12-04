% Do scenario MPC on the cart pendulum system
clear all

%% Setup
load('sysDyn.mat')
Q = [1 0 0 0;
     0 0 0 0;
     0 0 1 0;
     0 0 0 0];
R = 0.05;
sigma = 0.20; %[0.10 0.15 0.20 0.25];
M = length(sigma);
K = 4;
N = 10;
T = 150;
x0 = [1 0 0 0]'; % start the cart away from the origin

% allocate space to store closed loop path
xCL = NaN(4,T+1,M);
xOL = NaN(4,N+1,K,T);
xCL(:,1,:) = repmat(x0, 1, 1, M);
uCL = NaN(T,M);
uOL = NaN(N,T);

%% Solve
for m = 1:M
    for t = 1:T
        display([m t])
        omega = sigma(m)*randn(N,K);
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
        
        if ~strcmp(cvx_status, 'Solved')
            break
        end

        % store closed loop results
        uCL(t,m) = u(1);
        uOL(:,t) = u;
        xCL(:,t+1,m) = A*xCL(:,t,m) + B*(u(1) + sigma(m)*randn);
        xOL(:,:,:,t) = x;
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
% numViol = sum(abs(xCL(3,:))>0.1)
% save('resultsSigma2.mat')
