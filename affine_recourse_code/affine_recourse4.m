function [F,u0,status] = affine_recourse4(A,B,x_0,x_tar,n,T,mu,Sigma_bar,sig)
% Design of affine controller
% Adapted from example in "Design of Affine Controllers via Convex Optimization",
% Joelle Skaf and Stephen Boyd.
% www.stanford.edu/~boyd/papers/affine_contr.html
% Requires CVX package to run.
% @param A,B,x_0 -- gotten from the system dynamics, x_0 is initial state
% @param x_tar -- point we would like to track
% @param n = dimension of state space 
% @param T = Time horizon 0 - T;
% @param mu = mean of Gaussian noise, always zero for now
% @param Sigma_bar = covariance matrix (nxn)
% @param sig = standard deviation 

rand('state',0);
rand('state',0);

%Generate data matrices
%-----------------------
x0 = zeros(n*(T+1),1);
tmp = x_0;
for i=0:T
    x0(n*i+1:n*(i+1)) = tmp;
    tmp = A*tmp;
end

G = zeros(n*(T+1),n*T);
tmp = eye(n);
G(n+1:2*n,1:n) = tmp;
for i=2:T
    tmp = A*tmp;
    G(n*i+1:n*(i+1), 1:n*i) = [tmp G(n*(i-1)+1:n*i, 1:n*(i-1))];
end

H = zeros(n*(T+1),n*T);
tmp = B;
H(n+1:2*n,1:n) = tmp;
for i=2:T
    tmp = A*tmp;
    H(n*i+1:n*(i+1), 1:n*i) = [tmp H(n*(i-1)+1:n*i, 1:n*(i-1))];
end

S = zeros(n*T);
for i=0:T-1
    S(n*i+1:n*i+n, n*i+1:n*i+n) = Sigma_bar;
end


%Create affine controller
%------------------------
I = eye(n*(T+1));

cvx_begin quiet
variable Q(n*T, n*(T+1))
variable r(n*T)


Pxw = (I+H*Q)*G;
Puw = Q*G;
x_tilde = (I + H*Q)*x0 + H*r;
u_tilde = Q*x0 + r;


Z = Pxw*S;
U = Puw*S;
obj = 0;
for t=0:T
    
    Zt = Z(n*t+1:n*(t+1), :);
    obj = obj + square_pos(norm(Zt,'fro')) + square_pos(norm(x_tilde(n*t+1:n*(t+1), :) -x_tar)) ;
    
    if t ~= T
    Ut = U(n*t+1:n*(t+1), :);
    obj = obj + square_pos(norm(Ut,'fro')) + square_pos(norm(u_tilde(n*t+1:n*(t+1), :)));
    end

end

minimize (obj)

subject to

% constraints

for t=0:T-1
    
    
    P_x = Pxw(n*(t+1)+1:n*(t+2),:);
    P_u = Puw(n*t+1:n*(t+1),:);
    ut = u_tilde(n*t+1:n*(t+1), :);
    xt2 = x_tilde(n*(t+1)+1:n*(t+2), :);
    
    for i = 1:n
        
        -sig*norm(P_x(i,:),1) + xt2(i,:) >= 1;
        -sig*norm(P_u(i,:),1) + ut(i,:) >= -5;
        sig*norm(P_u(i,:),1) + ut(i,:) <= 5;
        
    end
    

end


% Q (n,n) block lower triangular
for i=0:T-1
    Q(n*i+1:n*(i+1), n*(i+1)+1:end) == 0
end

cvx_end

status = cvx_status;

I = eye(n*T);
F = (I+Q*H)\Q;
u0 = (I+Q*H)\r;
end

