function [F,u0,status] = affine_recourse4(A,B,x_0,x_tar,n,T,mu,Sigma_bar,sig)
% Dynamic portfolio optimization: income investing
% Example from "Design of Affine Controllers via Convex Optimization",
% Joelle Skaf and Stephen Boyd.
% www.stanford.edu/~boyd/papers/affine_contr.html
% Requires CVX package to run.
%clear all;


rand('state',0);
rand('state',0);

% data generation
% ---------------
% T = 10;                 % number of periods
% n = 2;                  % number of assets
% kappa = 0.1;            % transaction cost rate
% rho = .8;
% r_min = 0.10;

% return r(t) ~ exp(N(mu, Sigma))-1
% beta = 0.3+0.7*rand(n,1);           % beta's are unif on [0.3,1]
% beta(1) = 0;                        % beta for risk-free asset
% sigma_m = 0.1;                      % market std dev
% sigmas = 0.2*rand(n,1);             % firm-specific std dev unif on [0, 20]
% sigmas(1) = 0;                      % std dev for risk-free asset
% Sigma = sigma_m^2*beta*beta' + diag(sigmas.^2);

% mu_rf = 0.02;                       % risk-free return
% SR = 0.4;                           % sharpe ratio
% mu = mu_rf + SR*sqrt(diag(Sigma));  % expected return is 0.4*std dev

% [mu, idx] = sort(mu, 'ascend');
% Sigma = sigma_m^2*beta(idx)*beta(idx)' + diag(sigmas(idx).^2);

% rbar = exp(mu+0.5*diag(Sigma));
% Sigma_bar = (rbar*rbar').*(exp(Sigma)-1);
% rbar = rbar - 1;
% mu = [0;0];
% Sigma_bar = 0.1*eye(2);

% target portfolio
% ----------------
%xtar = [1;0;1;1]/3;
%x_tar = [1;1];

% affine controller
% -----------------
%I = eye(n);
%A = I + diag(rbar);
%B = I + diag(rbar);
%x_0 = xtar;
%A = [0.7 -0.1*2;-0.1*3 0.9];
%B = eye(2);
%x_0 = x_tar;

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
    %S(n*i+1:n*i+n, n*i+1:n*i+n) = diag(xtar)*sqrtm(Sigma_bar);
    S(n*i+1:n*i+n, n*i+1:n*i+n) = Sigma_bar;
end


% generate M random samples of r
%M = 1000;  % used in paper.
% M = 19;
% W = [];
% for t=0:T-1
%     %R = exp(mu*ones(1,M) + sqrtm(Sigma)*randn(n,M))-1;
%     %W = [W; (R - rbar*ones(1,M)).*(xtar*ones(1,M))];
%     R = mvnrnd(mu',Sigma_bar,M);
%     R = R';
%     W = [W;R];
%     
% end
%M = 1;

I = eye(n*(T+1));

cvx_begin
variable Q(n*T, n*(T+1))
variable r(n*T)
%variable W(n*T)
%variable x(n*(T+1))
%variable u(n*T)

Pxw = (I+H*Q)*G;
Puw = Q*G;
x_tilde = (I + H*Q)*x0 + H*r;
u_tilde = Q*x0 + r;

% x = Pxw*W + x_tilde*ones(1,M);
% u = Puw*W + u_tilde*ones(1,M);

Z = Pxw*S;
obj = 0;
for t=0:T
    
    Zt = Z(n*t+1:n*(t+1), :);
    obj = obj + square_pos(norm(Zt,'fro')) + square_pos(norm(x_tilde(n*t+1:n*(t+1), :) -x_tar));
    
%     P_x = Pxw(n*(t+1)+1:n*(t+2),:);
%     P_u = Puw(n*t+1:n*(t+1),:);
%     ut = u_tilde(n*t+1:n*(t+1), :);
%     xt2 = x_tilde(n*(t+1)+1:n*(t+2), :);
%     
%     temp_x = -inv((P_x + P_u))*(P_x*xt2 + P_u*ut);
%     
%     if(norm(temp_x,Inf) <= sig)
%         
%         temp_obj = -inv((P_x + P_u))*(P_x*xt2 + P_u*ut)+(xt2'*xt2)+(ut'*ut);
%         
%     else
%         w_low = -sig*ones(size(xt2));
%         w_high = sig*ones(size(xt2)); 
%         low = w_low'*(P_x + P_u)*w_low + 2*w_low'*(P_x*xt2 + P_u*ut) + (xt2'*xt2)+(ut'*ut);
%         high = w_high'*(P_x + P_u)*w_high + 2*w_high'*(P_x*xt2 + P_u*ut) + (xt2'*xt2)+(ut'*ut);
%      
%         temp_obj = max(low,high);
%     end
%         
%         
%     %obj = obj -square_pos (r_min+sum(ut) + kappa*norms(ut,1));
%     obj = obj + temp_obj;
%  

end

minimize (obj)

subject to

% constraints
%x - x_tilde <= Pxw*S*ones(n*T,1) ; %*ones(1,M);
%-x - x_tilde <= Pxw*S*ones(n*T,1);



%u <= Puw*(S.^(1/2))*ones(n*T,1) + u_tilde; %*ones(1,M);
%-u <= Puw*(S.^(1/2))*ones(n*T,1) + u_tilde;

for t=0:T-1
    
    %Zt = Z(n*t+1:n*(t+1), :);
    %square_pos(norm(Zt,'fro')) + square_pos(norm(x_tilde(n*t+1:n*(t+1), :) -xtar)) <= n*rho^2
    %sum(x_tilde(n*t+1:n*(t+1), :)) == 1;
    
    P_x = Pxw(n*(t+1)+1:n*(t+2),:);
    P_u = Puw(n*t+1:n*(t+1),:);
    ut = u_tilde(n*t+1:n*(t+1), :);
    xt2 = x_tilde(n*(t+1)+1:n*(t+2), :);
    
    for i = 1:n
        
        -sig*norm(P_x(i,:),1) + xt2(i,:) >= 1;
        -sig*norm(P_u(i,:),1) + ut(i,:) >= -5;
        sig*norm(P_u(i,:),1) + ut(i,:) <= 5;
        
    end
    
%     x >= repmat(x_tar*ones(1,M),T+1,1);
%     u <= repmat([5;5]*ones(1,M),T,1);l
%     -u <=repmat([5;5]*ones(1,M),T,1);

end

%W <= S*ones(n*T,1);
%-W <= S*ones(n*T,1);

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

