function c = solveTubeFHOCP(xk, zk, Phi, B, N, Nhat, lambda, V, g, h, bq)
% Solve the finite horizon optimal control problem 
% 
%
% xk: initial xcondition
% theta: sampled uncertainty for multiplicative disturbance (N x K)
% omega: sampled uncertainty for additive disturbance (2 x N x K)

% %% Problem data
% lambda = 0.422;
% V = [0.786  0.150;...   
%      0.150  0.072];
%  
% Nhat = 7;
% 
% z(:,1) = zk;
% 
% g = [1 0.2]';
% h = 1.2;
% 
% % bq = 0.086*ones(1,Nhat+N);
% % bq(1:7) = [0.013 0.055 0.073 0.080 0.083 0.085 0.085]; 
% 
% bq = 0.5*0.0585*ones(1,Nhat+N);
% bq(1:3) = 0.5*[0.0128 0.0535 0.0583];
% 

%% Solve problem
cvx_begin quiet
    variables c(1,N-1)
    variables z(2,N)
    Jk = [xk' c 1] * [xk' c 1]';
    minimize(Jk)
    subject to
        z(:,1) == zk;
        for i = 1:N-1
             z(:,i+1) == Phi*z(:,i) + B * c(:,i); %System Dynamics
             g'*z(:,i+1) <= h - bq(i+1)^0.5 * (g'*V^(-1)*g)^0.5;
        end
        
        for j = 1:Nhat
            g'*Phi^j <= h - bq(N+j)^0.5 * (g'*inv(V)*g)^0.5;
        end
cvx_end

end
