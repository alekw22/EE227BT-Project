function u = solveFHOCP_distOnlyExp(x0, omega1, omega2)
% Solve the finite horizon optimal control problem (FHOCP) in the numerical
% example from Schildbach et al (2014) with only the additive disturbance
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
cvx_begin
    variables u(2,N)
    for k = 1:K1
        x1_s1 = A*x0 + B*u(:,1) + omega1(:,1,k);
        x2_s1 = A*x1_s1 + B*u(:,2) + omega1(:,2,k);
        x3_s1 = A*x2_s1 + B*u(:,3) + omega1(:,3,k);
        x4_s1 = A*x3_s1 + B*u(:,4) + omega1(:,4,k);
        x5_s1 = A*x4_s1 + B*u(:,5) + omega1(:,5,k);
        objfun = objfun + x0'*x0 + u(:,1)'*u(:,1) + x1_s1'*x1_s1 + u(:,2)'*u(:,2) + ...
            x2_s1'*x2_s1 + u(:,3)'*u(:,3) + x3_s1'*x3_s1 + u(:,4)'*u(:,4) + x4_s1'*x4_s1 + ...
            u(:,5)'*u(:,5) + x5_s1'*x5_s1 + u(:,6)'*u(:,6);
    end
    for k = 1:K2
        x1_s2 = A*x0 + B*u(:,1) + omega2(:,1,k);
        x2_s2 = A*x1_s2 + B*u(:,2) + omega2(:,2,k);
        x3_s2 = A*x2_s2 + B*u(:,3) + omega2(:,3,k);
        x4_s2 = A*x3_s2 + B*u(:,4) + omega2(:,4,k);
        x5_s2 = A*x4_s2 + B*u(:,5) + omega2(:,5,k);
        objfun = objfun + x0'*x0 + u(:,1)'*u(:,1) + x1_s2'*x1_s2 + u(:,2)'*u(:,2) + ...
            x2_s2'*x2_s2 + u(:,3)'*u(:,3) + x3_s2'*x3_s2 + u(:,4)'*u(:,4) + x4_s2'*x4_s2 + ...
            u(:,5)'*u(:,5) + x5_s2'*x5_s2 + u(:,6)'*u(:,6);
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
