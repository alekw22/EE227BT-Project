% Runs the Tube MPC algorithm on the Schildbach (2014) example

clear; close all;
load sysDyn
load cart_pend

T2 = 150; 

Q = zeros(4);
Q(1,1) = 1;
Q(3,3) = 1;
K1 = 4;

R =0.05;

K = -dlqr(A,B,Q,R);
Phi = A + B*K;

W = 1;
sigma = 0.05;            %Noise

N = 6;                   %prediction horizon
Nhat = 7;

g1 = [0 0 1 0]';
h1 = 0.1;
lambda1 =  0.9366;
V1 = [ 136.0145   49.2399 -221.6827  -26.5812;...
       49.2399   54.8285 -117.9151  -33.3911;...
       -221.6827 -117.9151  765.9573   73.8023;...
     -26.5812  -33.3911   73.8023   29.8107];

bq1 = [0.0142 0.0252 0.0344 0.0426 0.0493 0.0554];

g2 = [0 0 -1 0]';
h2 = 0.1;
lambda2 = lambda1;
V2 = V1;
bq2 = bq1;

%Run 150 Steps
for k = 1:5
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
        u(:,i) = K*x(:,i) + c(:,1);
        
        w(:,i) = (.3)*rand(1) - 0.15;
        z(:,i+1) = Phi*z(:,i) + B*c(:,1);
        e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
        x(:,i+1) = z(:,i+1) + e(:,i+1);

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

toc