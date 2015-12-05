% Runs the Tube MPC algorithm on the cart pendulum example
clear; close all;
tic

T2 = 150; 
load sysDyn
load cartPendParam

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
n_star = 0;

g1 = [0 0 1 0]';
h1 = 0.1;
bq1 = [0.0142 0.0252 0.0344 0.0426 0.0493 0.0554];

g2 = [0 0 -1 0]';
h2 = 0.1;
bq2 = bq1;


for k = 1:10
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
    
    Cost(k,:) = C;
    cost(k) = sum(C);
    v1(k,:) = viol1;
    v2(k,:) = viol2;
end

toc