clear; close all;

T = 20;
T2 = 100;

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

N = 10;                 %prediction horizon
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


% bq1 = 0.01*bq1;           
% bq2 = 0.01*bq2;

%Initial Conditions
x = zeros(2,N);
e = zeros(2,N);
z = zeros(2,N);

x(:,1) = [1 1]';
e(:,1) = [0 0]';
z(:,1) = x(:,1);

u = zeros(2, N);

%Repeated Trials at k = 0
for j = 1:T
    s(j) = 0;
    %c = solveTubeFHOCP_S(x(:,1),z(:,1), Phi, B, N, Nhat, lambda1, V1, g1, h1, bq1);
    c = solveTubeFHOCP_S(x(:,1),z(:,1), Phi, B, N, Nhat, V1, g1, h1, bq1, V2, g2, h2, bq2);
    u(:,1) = K*x(:,1) + c(:,1);

    for i = 1:N-1;
        %w(:,i) = (0.224 + 0.224)*rand(2,1) - 0.224;
        w(:,i) = sigma*randn(2,1);
        z(:,i+1) = Phi*z(:,i) + B*c(:,i);
        e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
        x(:,i+1) = z(:,i+1) + e(:,i+1); 
        if (g1'*x(:,i+1) <= h1)
            s(j) = s(j)+1;
        end   
    end
    
figure(1)
hold on
plot(x(1,:)',x(2,:)')
end

%plot([-5 2],[(- g2(1)*(-5) + h2)/g2(2), (- g2(1)*(2) + h)/g(2)])

%Moving forward in time
x = zeros(2,T2-1);
e = zeros(2,T2-1);
z = zeros(2,T2-1);

x(:,1) = [1 1]';
e(:,1) = [0 0]';
z(:,1) = x(:,1);

u = zeros(2, T2);

for i = 1:T2
    c = solveTubeFHOCP_S(x(:,1),z(:,1), Phi, B, N, Nhat, V1, g1, h1, bq1, V2, g2, h2, bq2);
   % c = solveTubeFHOCP_S(x(:,i),z(:,i), Phi, B, N, Nhat,lambda1, V1, g1, h1, bq1);

    u(:,1) = K*x(:,1) + c(:,1);

    w(:,i) = (0.224 + 0.224)*rand(2,1) - 0.224;
    z(:,i+1) = Phi*z(:,i) + B*c(:,1);
    e(:,i+1) = Phi*e(:,i) + Bw*w(:,i);
    x(:,i+1) = z(:,i+1) + e(:,i+1); 
end

figure(2)
hold on
plot(x(1,:)',x(2,:)')

