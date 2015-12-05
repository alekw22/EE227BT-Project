function [F,u0,status] = affine_recourse_cart(A,B,Bw,x_0,x_tar,n,T,mu,Sigma_bar,sig,gamma)
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
nw = size(Bw,2);
G = zeros(n*(T+1),nw*T);
tmp = Bw;
G(n+1:2*n,1:nw) = tmp;
for i=2:T
    tmp = A*tmp;
    G(n*i+1:n*(i+1), 1:nw*i) = [tmp G(n*(i-1)+1:n*i, 1:nw*(i-1))];
end

nn = size(B,2);
H = zeros(n*(T+1),nn*T);
tmp = B;
H(n+1:2*n,1:nn) = tmp;
for i=2:T
    tmp = A*tmp;
    H(n*i+1:n*(i+1), 1:nn*i) = [tmp H(n*(i-1)+1:n*i, 1:nn*(i-1))];
end


S = zeros(nw*T);
for i=0:T-1
    S(nw*i+1:nw*i+nw, nw*i+1:nw*i+nw) = Sigma_bar;
end


%Create affine controller
%------------------------
I = eye(n*(T+1));

cvx_begin 
variable Q(nn*T, n*(T+1))
variable r(nn*T)


Pxw = (I+H*Q)*G;
Puw = Q*G;
x_tilde = (I + H*Q)*x0 + H*r;
u_tilde = Q*x0 + r;


Z = Pxw*0.15;
U = Puw*0.15;
obj = 0;
 R = eye(4,4);
 R(1,1) = 1; %1.2
 R(3,3) = 1; %1.1
 R(2,2) = 0; %0.3
 R(4,4) = 0; %).8
 R_u = 0.05;
for t=0:T
    
    Zt = Z(n*t+1:n*(t+1), :);
    obj = obj + square_pos(norm(Zt,'fro')) + square_pos(norm(R.^(1/2)*(x_tilde(n*t+1:n*(t+1), :) -x_tar))) ;
    
    
    if t ~= T
    Ut = U(nn*t+1:nn*(t+1), :);
    obj = obj + square_pos(norm(Ut,'fro')) + R_u*square_pos(norm(u_tilde(nn*t+1:nn*(t+1), :)));
    end
    
    obj = gamma*obj;

end


minimize (obj)

subject to

% constraints

g = [1;0.2];
for t=0:T-1
    
    
    P_x = Pxw(n*(t+1)+1:n*(t+2),:);
    P_u = Puw(nn*t+1:nn*(t+1),:);
    ut = u_tilde(nn*t+1:nn*(t+1), :);
    xt2 = x_tilde(n*(t+1)+1:n*(t+2), :);
    
    %sig*norm(g'*P_x,1) + g'*xt2 <= 1.2;
    
    
     for i = 3
%         
         -0.12*norm(P_x(i,:),1) + xt2(i,:) >= -0.1;
          0.12*norm(P_x(i,:),1) + xt2(i,:) <= 0.1;
%         %-sig*norm(P_u(i,:),1) + ut(i,:) >= -5;
%         %sig*norm(P_u(i,:),1) + ut(i,:) <= 5;
%         
    end
    

end


% Q (n,n) block lower triangular
for i=0:T-1
    Q(nn*i+1:nn*(i+1), n*(i+1)+1:end) == 0
end

cvx_end

status = cvx_status;

I = eye(nn*T);
F = (I+Q*H)\Q;
u0 = (I+Q*H)\r;
end

