function c = solveTubeFHOCP_S(xk, zk, Phi, B, N, Nhat, V1, g1, h1, bq1, V2, g2, h2, bq2)
%function c = solveTubeFHOCP(xk, zk, Phi, B, N, Nhat, lambda1, V1, g1, h1, bq1)

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
    variables c(2,N-1)
    variables z(2,N)
    Jk = [xk' c(1)' 1] * [xk' c(1)' 1]';
    minimize(Jk)
    subject to
        z(:,1) == zk;
        for i = 1:N-1
             z(:,i+1) == Phi*z(:,i) + B * c(:,i); %System Dynamics
             g1'*z(:,i+1) <= h1 - bq1(i+1)^0.5 * (g1'*V1^(-1)*g1)^0.5;
             g2'*z(:,i+1) <= h2 - bq2(i+1)^0.5 * (g2'*V2^(-1)*g2)^0.5;
        end
        
        for j = 1:Nhat
            g1'*Phi^j <= h1 - bq1(N+j)^0.5 * (g1'*inv(V1)*g1)^0.5;
            g2'*Phi^j <= h2 - bq2(N+j)^0.5 * (g2'*inv(V2)*g2)^0.5;
        end
        
        
cvx_end

end
