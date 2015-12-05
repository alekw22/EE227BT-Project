clear; close all;
tic
T1 = 200;

%System Model
A = [0.7 -0.2; -0.3 0.9];   
B = eye(2);
Q = eye(2);
R = eye(2);
K = -dlqr(A,B,Q,R);    %LQG optimal feedback
Phi = A + B*K;

Bw = eye(2);
W = eye(2);

sigma = 0.1;            %Noise

N = 6;                  %prediction horizon
Nhat = 7;  
n_star = 0;

g1 = [-1 0]';
h1 = -1;
lambda1 =  0.3234;
V1 = calcTubeParamV(lambda1);
bq1 = 0.8859*ones(1,Nhat+N);
bq1(1:4) =[ 0.3219 0.8380 0.8380 0.8859];

g2 = [0 -1]';
h2 = -1;
lambda2 = 0.3728;
V2 = calcTubeParamV(lambda2);
bq2 = 0.9815*ones(1,Nhat+N);
bq2(1:4) = [0.3214 0.9098 0.9614 0.9815];

%Run 6 Steps 200 times

for k = 1:200
    k
    x = zeros(2,N);
    e = zeros(2,N);
    z = zeros(2,N);

    x(:,1) = [1 1]';
    e(:,1) = [0 0]';
    z(:,1) = x(:,1);

    u = zeros(2, N);

    for i = 1:N
        c = solveTubeFHOCP_S(x(:,i),z(:,i), Phi, B, N, Nhat, V1, g1, h1, bq1, V2, g2, h2, bq2);
        u(:,i) = K*e(:,i) + c(:,1);

        w(:,i) = sigma^0.5*randn(2,1);
        z(:,i+1) = Phi*z(:,i) + B*c(:,1);
        e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
        %x(:,i+1) = z(:,i+1) + e(:,i+1); 
        
        C(i) = norm(u(:,i))^2 + norm(x(:,i))^2;
        viol1(i) = (g1'*x(:,i) > h1);
        viol2(i) = (g2'*x(:,i) > h2); 
    end
    figure(1)
    hold on
    plot(x(1,:)',x(2,:)')
    Cost(k,:) = C;
    cost(k) = sum(C);
    v1(k,:) = viol1;
    v2(k,:) = viol2;
end


%Run 3000 Steps

% for k = 1:3000
%     
%     x = zeros(2,T2);
%     e = zeros(2,T2);
%     z = zeros(2,T2);
% 
%     x(:,1) = [1 1]';
%     e(:,1) = [0 0]';
%     z(:,1) = x(:,1);
% 
%     u = zeros(2,T2);
% 
%     for i = 1:T2
%         i
%         c = solveTubeFHOCP_S(x(:,i),z(:,i), Phi, B, N, Nhat, V1, g1, h1, bq1, V2, g2, h2, bq2);
%         u(:,i) = K*x(:,i) + c(:,1);
% 
%         w(:,i) = sigma^0.5*randn(2,1);
%         z(:,i+1) = Phi*z(:,i) + B*c(:,1);
%         e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
%         x(:,i+1) = z(:,i+1) + e(:,i+1); 
%         
%         C(i) = norm(u(:,i))^2 + norm(x(:,i))^2;
%         viol1(i) = (g1'*x(:,i) > h1);
%         viol2(i) = (g2'*x(:,i) > h2); 
%     end
%     figure(1)
%     hold on
%     plot(x(1,:)',x(2,:)')
%     Cost(k,:) = C;
%     cost(k) = sum(C);
%     v1(k,:) = viol1;
%     v2(k,:) = viol2;
% end
% 
% figure(2)
% hold on
% plot(x(1,:)',x(2,:)')
% xlim([0,4])
% ylim([0,4])
 
toc