function c = solveTubeFHOCP(xk, zk, Phi, B, N, Nhat, lambda, V, g, h, bq)
% Solve the finite horizon optimal control problem 
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
