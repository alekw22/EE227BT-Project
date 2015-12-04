% Tube scaling distributions for cart-pendulum system
clear

%% problem data
wmax = 0.15;
alpha_bar = wmax^2;
lambda = 0.9366;                   
beta_bar = alpha_bar/(1-lambda);     % maximum error tube scaling
rho = 1000;                          % number of discrete samples
h = beta_bar/rho;
omega = ceil(beta_bar*(1-lambda)/(lambda*h));    % extra discrete samples to compute CDF

%% P matrix for advancing the distribution
P = zeros(rho+1, rho+omega+1);
x = linspace(0, beta_bar/lambda, rho+omega+1);
for i = 1:(rho+1)
    i
    for j = 1:(rho+omega)
        midpoint = x(j) + h/2;
        P(i,j) = lambda*h*f_alpha(x(i) - lambda*midpoint, wmax);
    end
end

advanceCDF = @(pi,mu) P*[pi; ones(omega,1)];
mu = 100;

%%
n = 15;
x = linspace(0,beta_bar,rho+1);
piU = cell(n);         % cell array of CDFs
piS = cell(n);         % cell array of smoothed CDFs
pi0 = ones(rho+1, 1);   
piU{1,1} = advanceCDF(pi0);
piS{1,1} = smoothCDF(piU{1,1}, mu);

for k = 1:(n-1)
    k
    piU{k+1,1} = advanceCDF(piS{k,1});
    piS{k+1,1} = smoothCDF(piU{k+1,1}, mu);
end

%%
figure
plot(x, pi0)
hold on
for k = 1:n
    plot(x, piS{k,1})
end

%% find b_bar
p = 0.8;
q = 2*p - 1;

b = zeros(n,1);

for k = 1:n
    b(k) = x(find(piS{k,1} > q, 1));
end
