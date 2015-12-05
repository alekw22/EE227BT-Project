function cvx_optval = calcTubeParam(lambda)

%Cannon
% g = [1 0.2]';
% A = [ 1.6   1.1; ...
%      -0.7   1.2];
% B = [1 1]';
% Q = eye(2);
% R = 1;
% K1=2;

%Schil
% A = [0.7 -0.2; -0.3 0.9];   
% B = eye(2);
% %g1 = [-1 0]';
% g = [0 -1]';
% Q = eye(2);
% R = eye(2);
% K1 = 2;

%Inverted Pendulum
load sysDyn

Q = zeros(4);
Q(1,1) = 1;
Q(3,3) = 1;
K1 = 4;

R =0.05;
g = [0 0 1 0]';
%g = [0 0 -1 0]';

K = -dlqr(A,B,Q,R);
Phi = A + B*K;

%Bw = eye(2);
W = 1;

cvx_begin quiet
variables U(K1,K1) 
minimize((1/(1 - lambda)*g'*U*g))
subject to:
    (U - 1/lambda*Phi*U*Phi' - Bw*W*Bw') == semidefinite(K1)
    U == semidefinite(K1)
cvx_end

V = inv(U);

end