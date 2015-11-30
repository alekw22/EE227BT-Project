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
T = 10;                 % number of periods
n = 2;                  % number of assets
kappa = 0.1;            % transaction cost rate
rho = .8;
r_min = 0.10;

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
mu = [0;0];
Sigma_bar = 0.1*eye(2);

% target portfolio
% ----------------
%xtar = [1;0;1;1]/3;
xtar = [1;1];

% affine controller
% -----------------
%I = eye(n);
%A = I + diag(rbar);
%B = I + diag(rbar);
%x_0 = xtar;
A = [0.7 -0.1*2;-0.1*3 0.9];
B = eye(2);
x_0 = xtar;

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
M = 100;
W = [];
for t=0:T-1
    %R = exp(mu*ones(1,M) + sqrtm(Sigma)*randn(n,M))-1;
    %W = [W; (R - rbar*ones(1,M)).*(xtar*ones(1,M))];
    R = mvnrnd(mu',Sigma_bar,M);
    R = R';
    W = [W;R];
    
end
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

x = Pxw*W + x_tilde*ones(1,M);
u = Puw*W + u_tilde*ones(1,M);


obj = 0;
for t=0:T-1
    ut = u(n*t+1:n*(t+1), :);
    xt2 = x(n*t+1:n*(t+1), :);
    %obj = obj -square_pos (r_min+sum(ut) + kappa*norms(ut,1));
    obj = obj + square_pos(norms(ut,2)) + square_pos(norms(xt2,2));
    
end
minimize (sum(obj)/(M*T))

subject to

% constraints
%x - x_tilde <= Pxw*S*ones(n*T,1) ; %*ones(1,M);
%-x - x_tilde <= Pxw*S*ones(n*T,1);



%u <= Puw*(S.^(1/2))*ones(n*T,1) + u_tilde; %*ones(1,M);
%-u <= Puw*(S.^(1/2))*ones(n*T,1) + u_tilde;
%Z = Pxw*S;
for t=0:T
    %Zt = Z(n*t+1:n*(t+1), :);
    %square_pos(norm(Zt,'fro')) + square_pos(norm(x_tilde(n*t+1:n*(t+1), :) -xtar)) <= n*rho^2
    
    %sum(x_tilde(n*t+1:n*(t+1), :)) == 1;
    x >= repmat(xtar*ones(1,M),T+1,1);
    u <= 2*repmat([5;5]*ones(1,M),T,1);
    -u <= 2*repmat([5;5]*ones(1,M),T,1);
end

%W <= S*ones(n*T,1);
%-W <= S*ones(n*T,1);

% Q (n,n) block lower triangular
for i=0:T-1
    Q(n*i+1:n*(i+1), n*(i+1)+1:end) == 0
end
cvx_end

I = eye(n*T);
F = (I+Q*H)\Q;
u0 = (I+Q*H)\r;
%%
%test affine controller
%----------------------
%generate m random samples of r
m = 1;
Wt = [];
Rs = [];
for t=0:T-1
    %     R = exp(mu*ones(1,m) + sqrtm(Sigma)*randn(n,m))-1;
    %     Rs = [Rs; R];
    %     Wt = [Wt; (R - rbar*ones(1,m)).*(xtar*ones(1,m))];
    
    R = mvnrnd(zeros(size(mu')),zeros(size(Sigma_bar)),m);
    R = R';
    Wt = [Wt;R];
end

% propagating using linearized dynamics
x_linear = Pxw*Wt + x_tilde*ones(1,m);
u_linear  = Puw*Wt + u_tilde*ones(1,m);
%util_linear = 0;
%cash_linear = zeros(T,m);
% for t=0:T-1
%     ut = u_linear(n*t+1:n*(t+1), :);
%     util_linear = util_linear -square_pos (r_min+sum(ut) + kappa*norms(ut,1));
%     cash_linear(t+1,:) = -sum(ut) - kappa*norms(ut,1);
% end
% util_linear = util_linear/T;
% cash_linear = -sum(u_linear) - kappa*norms(u_linear,1);

%%
%plot result
X_new_1 = x_linear(1:2:end,:);
X_new_2 = x_linear(2:2:end,:);

%plot(X_new_1',X_new_2');

for j = 1:size(X_new_1,2)
    
    plot(X_new_1(:,j), X_new_2(:,j));
    hold on;
    
    
end
%%

x1 = get(gca,'xlim');
plot(x1,[1,1]);

y1 = get(gca,'ylim');
plot([1,1],y1);

