clear

%% Problem data
A = [0.7 -0.25; -0.4 0.9];
B = eye(2);

N = 5;       % prediction horizon
x0 = [2; 2]; % initial condition

%% Solve problem
cvx_begin
    variables x(2,N+1) u(2,N)
    objfun = 0;
    for k = 1:5
        objfun = objfun + x(:,k)'*x(:,k) + u(:,k)'*u(:,k);
    end
    minimize( objfun )
    subject to
        -5 <= u <= 5   % input constraints
        1 <= x         % state constraints
        x(:,1) == x0   % initial condition constraint
        for k = 1:N
            x(:,k+1) == A*x(:,k) + B*u(:,k)  % system dynamics constraints
        end
cvx_end

%% Plot output
figure
plot(x(1,:), x(2,:))
xlabel('x1'); ylabel('x2'); title('Open loop trajectory')
