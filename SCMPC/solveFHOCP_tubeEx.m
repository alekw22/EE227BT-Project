function [u, x] = solveFHOCP_tubeEx(x0, omega)
% Solve the finite horizon optimal control problem (FHOCP) in the numerical
% example from Cannon et al (2010). 
%
% x0: initial condition
% omega: sampled uncertainty for disturbance (2 x N x K)

%% Problem data
A = [1.6 1.1; -0.7 1.2];
B = [1 1]';
g = [1 0.2]';
h = 1.2;
N = size(omega,2);    % prediction horizon
K = size(omega,3);    % number of scenarios  

%% Solve problem
cvx_begin quiet
    variables x(2,N+1,K) u(1,N)
    objfun = 0;
    for k = 1:K
        for i = 1:N
            objfun = objfun + x(:,i,k)'*x(:,i,k) + u(:,i)'*u(:,i);
        end
    end
    minimize(objfun)
    subject to
        for k = 1:K
            x(:,1,k) == x0;   % initial condition constraint
            for i = 1:N
                x(:,i+1,k) == A*x(:,i,k) + B*u(:,i) + omega(:,i,k); % system dynamics constraints
                g'*x(:,i+1,k) <= h;     % state constraint
            end
        end
cvx_end
