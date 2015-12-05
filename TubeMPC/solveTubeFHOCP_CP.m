function c = solveTubeFHOCP_CP(xk, zk, Phi, B, N, Nhat, V1, g1, h1, bq1, V2, g2, h2, bq2, Q, R, K)
% Solve the finite horizon optimal control problem 
% 
% xk: initial xcondition
% 

%% Solve problem
cvx_begin quiet
    variables c(1,N-1)
    variables z(length(xk),N)
   %Jk = [(xk - [-1 0 0 0 ]')' c(1)' 1] *[Q, zeros(4,2); zeros(1,4) R 0; zeros(1, 5) 1] * [(xk - [-1 0 0 0 ]')' c(1)' 1]';
   Jk = [(xk') c(1)' 1] *[Q, zeros(4,2); zeros(1,4) R 0; zeros(1, 5) 1] * [(xk)' c(1)' 1]';
   %Jk = 20000*(xk - [-1 0 0 0 ]')'*Q*(xk - [-1 0 0 0 ]') + (c(:,1))*R*(c(:,1));
   %Jk = norm(c); 
   minimize(Jk)
    subject to
        z(:,1) == zk;
        for i = 1:N-1
             z(:,i+1) == Phi*z(:,i) + B * c(:,i); %System Dynamics
             g1'*z(:,i+1) <= h1 - bq1(i+1)^0.5 * (g1'*V1^(-1)*g1)^0.5;
             g2'*z(:,i+1) <= h2 - bq2(i+1)^0.5 * (g2'*V2^(-1)*g2)^0.5;
        end
%         
%         for j = 1:Nhat
%             g1'*Phi^j <= h1 - bq1(N+j)^0.5 * (g1'*inv(V1)*g1)^0.5;
%             g2'*Phi^j <= h2 - bq2(N+j)^0.5 * (g2'*inv(V2)*g2)^0.5;
%         end
%         
        
cvx_end

end
