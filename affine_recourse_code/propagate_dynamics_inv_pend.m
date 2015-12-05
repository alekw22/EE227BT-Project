function [ res1,res2,res3 ] = propagate_dynamics_inv_pend(A,B,Bw,T,u,x_01,x_02,F,u_0,mu,Sigma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n = size(x_01,1);
nn = size(B,2);
nw = size(Bw,2);
x = zeros(n*(T+1),1); 
x_nfb = zeros(n*(T+1),1); 
U_fb = zeros(nn*T,1);
x(1:n,:) = x_01;
x_nfb(1:n,:) = x_02;

%generate m random samples of r
m = 1;
Wt = [];
for t=0:T-1
    %     R = exp(mu*ones(1,m) + sqrtm(Sigma)*randn(n,m))-1;
    %     Rs = [Rs; R];
    %     Wt = [Wt; (R - rbar*ones(1,m)).*(xtar*ones(1,m))];
    
    R = mvnrnd(mu,Sigma,m);
    R = R';
    Wt = [Wt;R];
end

for i = 1:T
    u_fb = u_0(nn*(i-1)+1:nn*i,:) + F(nn*(i-1)+1:nn*i,1:n*i)*x(1:n*i,:);
    x((n*i)+1:n*(i+1),:) = A*x(n*(i-1)+1:n*i,:) + B*u_fb + Bw*Wt(nw*(i-1)+1:nw*i,:);
    U_fb(nn*(i-1)+1:nn*i,:) = u_fb;
    
end

for i = 1:T
    %u_fb = u_0(n*(i-1)+1:n*i,:) + F(n*(i-1)+1:n*i,1:n*i)*x(1:n*i,:);
    x_nfb((n*i)+1:n*(i+1),:) = A*x_nfb(n*(i-1)+1:n*i,:) + B*u(nn*(i-1)+1:nn*i,:)+ Bw*Wt(nw*(i-1)+1:nw*i,:);  
    
end

% %plot result
% X_new_1 = x(1:2:end,:);
% X_new_2 = x(2:2:end,:);
% 
% X_nfb_1 = x_nfb(1:2:end,:);
% X_nfb_2 = x_nfb(2:2:end,:);
% 
% %plot(X_new_1',X_new_2');
% 
% for j = 1:size(X_new_1,2)
%     
%     plot(X_new_1(:,j), X_new_2(:,j));
%     hold on;
%       
% end
% hold off
% 
% figure
% for j = 1:size(X_nfb_1,2)
%     
%     plot(X_nfb_1(:,j), X_nfb_2(:,j));
%     hold on;
%       
% end
% title('no feedback');

res1 = x;
res2 = x_nfb;
res3 = U_fb;

end



