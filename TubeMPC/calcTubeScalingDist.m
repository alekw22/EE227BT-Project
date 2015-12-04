% Calculate the tube scaling distributions for the example of Cannon et al.
% (2010)
clear

%% problem data
alpha_bar = 0.1;         % maximum disturbance tube scaling
sigma = 1/12;            % standard deviation of disturbance
lambda = 0.422;                   
beta_bar = alpha_bar/(1-lambda);     % maximum error tube scaling
rho = 1000;                          % number of discrete samples
h = beta_bar/rho;
omega = ceil(beta_bar*(1-lambda)/(lambda*h));    % extra discrete samples to compute CDF

%% density of disturbance scaling
tail = 1 - chi2cdf(alpha_bar/(sigma^2),2);
scale = 1/(1-tail);
f_alpha = @(x) (x <= alpha_bar).*chi2pdf(x/(sigma^2),2)*scale/(sigma.^2);
F_alpha = @(x) min([scale*chi2cdf(x/(sigma^2),2); ones(1,length(x))]);

%% P matrix for advancing the distribution
load('matrixP_padded.mat')
% P = zeros(rho+1, rho+omega+1);
% x = linspace(0, beta_bar/lambda, rho+omega+1);
% for i = 1:rho
%     i
%     for j = 1:(rho+omega+1)
%         if j > 1 && j < (rho+omega+1)
%             P(i,j) = lambda*h*f_alpha(x(i) - lambda*x(j));
%         else
%             P(i,j) = lambda*(h/2)*f_alpha(x(i) - lambda*x(j));
%         end
%     end
% end
% P(end,end) = 1;

advanceCDF = @(pi,mu) P*[pi; ones(omega,1)];
mu = 100;

%% Calculate the distribution approximations given beta_0 = 0
n = 8;
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

%% Plot all pi_i|0
figure
plot(x, pi0)
hold on
for k = 1:n
    plot(x, piS{k,1})
end

%% Find b_bar, the maximum tube scalings at each time step
p = 0.8;
q = 2*p - 1;

b = zeros(n);

for k = 1:n
    b(k,1) = x(find(piS{k,1} > q, 1));
end

for i = 2:n
    if i > 2
        for j = 2:(i-1)
            disp([i j])
            piU{i,j} = advanceCDF(piS{i-1,j});
            piS{i,j} = smoothCDF(piU{i,j}, mu);
            b(i,j) = x(find(piS{i,j} > q, 1));
        end
    end
    
    indDet = find(piS{i-1,1} >= 0.9999, 1);
    piU{i,i} = P*[zeros(indDet,1); ones(rho+omega+1-indDet,1)];
    piS{i,i} = smoothCDF(piU{i,i}, mu);
    b(i,i) = x(find(piS{i,i} > q, 1));
end

b_bar = max(b,[],2);

