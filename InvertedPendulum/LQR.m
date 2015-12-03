% LQR on the inverted pendulum system
clear

%% Setup
load('sysDyn.mat')
Q = [1 0 0 0;
     0 0 0 0;
     0 0 1 0;
     0 0 0 0];
R = 1;
[K,S,E] = dlqr(A,B,Q,R);
Acl = A-B*K;

%% Simulate
T = 500;
x0 = [1 0 0 0]'; % start the cart away from the origin

% allocate space to store closed loop path
xCL = zeros(4,T+1);
xCL(:,1) = x0;
uCL = zeros(1,T);

for t = 1:T
    uCL(t) = -K*xCL(:,t);
    xCL(:,t+1) = Acl*xCL(:,t);
end

%% Plot
figure
subplot(211)
plot((0:T)*0.02, xCL(1,:))
ylabel('x')

subplot(212)
plot((0:T)*0.02, xCL(3,:))
ylabel('\phi')
xlabel('Time')

figure
plot((0:(T-1))*0.02, uCL)
ylabel('Input')
xlabel('Time')


