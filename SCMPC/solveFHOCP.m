function [u, x] = solveFHOCP(x0, theta, omega)
% Solve the finite horizon optimal control problem (FHOCP) in the numerical
% example from Schildbach et al (2014). 
%
% x0: initial condition
% theta: sampled uncertainty for multiplicative disturbance (N x K)
% omega: sampled uncertainty for additive disturbance (2 x N x K)

%% Problem data
Anom = [0.7 -0.2; -0.3 0.9];   % nominal state transition matrix
Aunc = [0 -0.1; -0.2 0];       % coefficients of random variable theta
B = eye(2);
N = size(theta,1);    % prediction horizon
K = size(theta,2);    % number of scenarios  

%% Solve problem
cvx_begin
    variables x(2,N+1,K) u(2,N)
    objfun = 0;
    for k = 1:K
        for i = 1:N
            objfun = objfun + x(:,i,k)'*x(:,i,k) + u(:,i)'*u(:,i);
        end
    end
    minimize(objfun)
    subject to
        -5 <= u <= 5          % input constraints
        1 <= x(:,2:end,:)         % state constraints
        for k = 1:K
            x(:,1,k) == x0;   % initial condition constraint
            for i = 1:N
                % system dynamics constraints
                x(:,i+1,k) == (Anom + Aunc*theta(i,k))*x(:,i,k) + B*u(:,i) ...
                    + omega(:,i,k)
            end
        end
cvx_end
