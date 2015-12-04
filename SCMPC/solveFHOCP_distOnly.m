function u = solveFHOCP_distOnly(x0, omega1, omega2)
% Solve the finite horizon optimal control problem (FHOCP) in the numerical
% example from Schildbach et al (2014) with only the additive disturbance.
% Treat the two state constraints individually.
%
% x0: initial condition
% omega: sampled uncertainty for additive disturbance (2 x N x K1+K2)
% K1: number of samples for constraint on x1
% K2: number of samples for constraint on x2

%% Problem data
A = [0.7 -0.2; -0.3 0.9];
B = eye(2);
N = size(omega1,2);    % prediction horizon
K1 = size(omega1,3);
K2 = size(omega2,3);

%% Solve problem
cvx_begin quiet
    variables x_s1(2,N+1,K1) x_s2(2,N+1,K2) u(2,N)
    objfun = 0;
    for k = 1:K1
        for i = 1:N
            objfun = objfun + x_s1(:,i,k)'*x_s1(:,i,k) + u(:,i)'*u(:,i);
        end
    end
    for k = 1:K2
        for i = 1:N
            objfun = objfun + x_s2(:,i,k)'*x_s2(:,i,k) + u(:,i)'*u(:,i);
        end
    end
    minimize(objfun)
    subject to
        -5 <= u <= 5              % input constraints
        1 <= x_s1(1,2:end,:)         % state constraint on x1
        1 <= x_s2(2,2:end,:)         % state constraint on x2
        for k = 1:K1
            x_s1(:,1,k) == x0;   % initial condition constraint
            for i = 1:N
                % system dynamics constraints
                x_s1(:,i+1,k) == A*x_s1(:,i,k) + B*u(:,i) + omega1(:,i,k)
            end
        end
        for k = 1:K2
            x_s2(:,1,k) == x0;   % initial condition constraint
            for i = 1:N
                % system dynamics constraints
                x_s2(:,i+1,k) == A*x_s2(:,i,k) + B*u(:,i) + omega2(:,i,k)
            end
        end
cvx_end
