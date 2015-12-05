% Runs the Tube MPC algorithm on the cart pendulum example on Cannon (2014)
clear; close all;

A = [ 1.6   1.1;  -0.7   1.2];  %System Model
B = [1 1]';

Q = eye(2);
R = 1;
K = -dlqr(A,B,Q,R);             %LQG optimal feedback

Phi = A + B*K;

Bw = eye(2);
W = eye(2);

g = [1 0.2]';
h = 1.2;
p = 0.1;
N = 6;                       %prediction horizon
Nhat = 7;  
n_star = 0;

lambda = 0.422;
V = [0.786  0.150;...   
     0.150  0.072];
 
bq = 0.086*ones(1,Nhat+N);
bq(1:7) = [0.013 0.055 0.073 0.080 0.083 0.085 0.085];  

%
% bq =[0.0128 0.0196 0.0221 0.0232 0.0235 0.0239];

for k = 1:200
    k
    x = zeros(2,N);
    e = zeros(2,N);
    z = zeros(2,N);

    x(:,1) = [-5 60]';
    e(:,1) = [0 0]';
    z(:,1) = x(:,1);

    u = zeros(2, N);

    for i = 1:N
        c = solveTubeFHOCP(x(:,i),z(:,i), Phi, B, N, Nhat, lambda, V, g, h, bq);
        u(1) = K*x(:,1) + c(1);
        
        %w(:,i) = (0.224 + 0.224)*rand(2,1) - 0.224;
        w(:,i) = 1/12*randn(2,1);
        z(:,i+1) = Phi*z(:,i) + B*c(:,1);
        e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
        x(:,i+1) = z(:,i+1) + e(:,i+1); 
        
        C(i) = norm(u(:,i))^2 + norm(x(:,i))^2;
        viol(i) = (g'*x(:,i) > h);
        
    end
    X{k} = x;
    figure(1)
    subplot(1,2,1)
    hold on
    plot(x(1,:)',x(2,:)')
    plot([-5 2],[(- g(1)*(-5) + h)/g(2), (- g(1)*(2) + h)/g(2)],'black')
    xlim([-5,2])
    ylim([-10,60])
    
    subplot(1,2,2)
    hold on
    plot(x(1,:)',x(2,:)','.', 'MarkerSize', 15)
    plot([-5 2],[(- g(1)*(-5) + h)/g(2), (- g(1)*(2) + h)/g(2)],'black')
    xlim([-2.5,-1.5])
    ylim([13,19])
    
    Cost(k,:) = C;
    cost(k) = sum(C);
    v(k,:) = viol;
end

% %Repeated Trials at k = 0
% for j = 1:T
%     s(j) = 0;
%     c = solveTubeFHOCP(x(:,1),z(:,1), Phi, B, N, Nhat, lambda, V, g, h, bq);
%     u(1) = K*x(:,1) + c(1);
% 
%     %noise
%     for i = 1:N-1;
%         %w(:,i) = (0.224 + 0.224)*rand(2,1) - 0.224;
%         w(:,i) = randn;
%         z(:,i+1) = Phi*z(:,i) + B*c(i);
%         e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
%         x(:,i+1) = z(:,i+1) + e(:,i+1); 
%         if (g'*x(:,i+1) <= h)
%             s(j) = s(j)+1;
%         end   
%     end
%     
% figure(1)
% hold on
% plot(x(1,:)',x(2,:)')
% end

plot([-5 2],[(- g(1)*(-5) + h)/g(2), (- g(1)*(2) + h)/g(2)])

% %Moving forward in time
% x = zeros(2,T-1);
% e = zeros(2,T-1);
% z = zeros(2,T-1);
% 
% x(:,1) = [-5 60]';
% e(:,1) = [0 0]';
% z(:,1) = x(:,1);
% 
% u = zeros(1, T);
% 
% for i = 1:T
%     c = solveTubeFHOCP(x(:,i),z(:,i), Phi, B, N, Nhat,lambda, V, g, h, bq);
% 
%     u(i) = K*x(:,i) + c(1);
% 
%     w(:,i) = (0.224 + 0.224)*rand(2,1) - 0.224;
%     z(:,i+1) = Phi*z(:,i) + B*c(1);
%     e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
%     x(:,i+1) = z(:,i+1) + e(:,i+1); 
% end
% 
% figure(2)
% hold on
% plot(x(1,:)',x(2,:)')
% 
