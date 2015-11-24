% Replicate the results of Cannon et al. (2010), more specifically the tube
% scaling distributions
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

advanceCDF = @(pi) smoothCDF(P*[pi; ones(omega,1)]);

%%
n = 8;
x = linspace(0,beta_bar,rho+1);
pi = cell(n);
pi0 = ones(rho+1, 1);   % cell array of CDFs
pi{1,1} = advanceCDF(pi0);

for k = 1:(n-1)
    k
    pi{k+1,1} = advanceCDF(pi{k,1});
end

%%
figure
plot(x, pi0)
hold on
for k = 1:n
    plot(x, pi{k,1})
end

%% find b_bar
p = 0.8;
q = 2*p - 1;

b = zeros(n);

for k = 1:n
    b(k,1) = x(find(pi{k,1} > q, 1));
end

for i = 2:n
    if i > 2
        for j = 2:(i-1)
            disp([i j])
            pi{i,j} = advanceCDF(pi{i-1,j});
            b(i,j) = x(find(pi{i,j} > q, 1));
        end
    end
    
    indDet = find(pi{i-1,1} >= 0.999999, 1);
    pi{i,i} = P*[zeros(indDet,1); ones(rho+omega+1-indDet,1)];
    pi{i,i} = smoothCDF(pi{i,i});
    b(i,i) = x(find(pi{i,i} > q, 1));
end

b_bar = max(b,[],2);

