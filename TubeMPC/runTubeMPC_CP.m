clear; close all;
tic
T1 = 200;
T2 = 150; 
%System Model
%Inverted Pendulum
load sysDyn
load cart_pend

Q = zeros(4);
Q(1,1) = 1;
Q(3,3) = 1;
K1 = 4;

R =0.05;

K = -dlqr(A,B,Q,R);
Phi = A + B*K;

W = 1;

sigma = 0.05;            %Noise

N = 6;                  %prediction horizon
Nhat = 7;  
n_star = 0;

g1 = [0 0 1 0]';
h1 = 0.1;
% lambda1 =  0.9366;
% V1 = calcTubeParamV(lambda1);
bq1 = 0.5331*ones(1,Nhat+N);
bq1(1:3) =[0.1167 0.4874 0.5315];
bq1 = 0.1*bq1;

bq1 = [0.0142 0.0252 0.0344 0.0426 0.0493 0.0554];

g2 = [0 0 -1 0]';
h2 = 0.1;
% lambda2 = 0.3728;
% V2 = calcTubeParamV(lambda2);
bq2 = bq1;

%Run 6 Steps 200 times

% for k = 1:200
%     k
%     x = zeros(2,N);
%     e = zeros(2,N);
%     z = zeros(2,N);
% 
%     x(:,1) = [1 1]';
%     e(:,1) = [0 0]';
%     z(:,1) = x(:,1);
% 
%     u = zeros(2, N);
% 
%     for i = 1:N
%         c = solveTubeFHOCP_S(x(:,i),z(:,i), Phi, B, N, Nhat, V1, g1, h1, bq1, V2, g2, h2, bq2);
%         u(:,i) = K*e(:,i) + c(:,1);
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


% Run 3000 Steps

coeff = [1 0.1 0.01];

for k = 1:1
    k
    x = zeros(4,T2);
    e = zeros(4,T2);
    z = zeros(4,T2);

    x(:,1) = [1 0 0 0]';
    e(:,1) = [0 0 0 0]';
    z(:,1) = x(:,1);

    u = zeros(1,T2);

    for i = 1:T2
        
        c = solveTubeFHOCP_CP(x(:,i),z(:,i), Phi, B, N, Nhat, V1, g1, h1, bq1, V2, g2, h2, bq2, Q, R, K);
        control(:,i) = c(:,1);

        u(:,i) = K*x(:,i) + c(:,1);
        
        w(:,i) = (.3)*rand(1) - 0.15;
        z(:,i+1) = Phi*z(:,i) + B*c(:,1);
       
        e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
        x(:,i+1) = z(:,i+1) + e(:,i+1);
       % x(:,i+1) = A*x(:,i) + B*u(:,i) + w(:,i);

        C(i) =u(:,i)*R*u(:,i) + x(:,i)'*Q*(x(:,i));
        
        viol1(i) = (g1'*x(:,i) > h1);
        viol2(i) = (g2'*x(:,i) > h2); 
    end
    
    X{k} = x;
    U{k} = u;
    
    figure(1)
    hold on
    plot(x(1,:)')
    plot(x(3,:)')
    Cost(k,:) = C;
    cost(k) = sum(C);
    v1(k,:) = viol1;
    v2(k,:) = viol2;
    
    
end

figure(5)
hold on
plot(x(1,:)')
plot(x(3,:)')
% xlim([0,4])
% ylim([0,4])
 
toc