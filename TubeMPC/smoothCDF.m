function CDFs = smoothCDF(CDF, mu)
% smooth a discretized approximation of a Cumulative Discretization
% Function subject to constraints that make it a CDF

rho = length(CDF) - 1;
Q = 2*eye(rho+1);
Q(1,1) = 1;
Q(end,end) = 1;
for k = 1:rho
    Q(k, k+1) = -1;
    Q(k+1, k) = -1;
end

cvx_begin quiet
    variable CDFs(rho+1)
    objfun = norm(CDFs - CDF) + mu*CDFs'*Q*CDFs;
    minimize(objfun)     % curve fit and promote smoothness
    subject to
    0 <= CDFs <= 1                 % must be a probability
    CDFs(1) == CDF(1);
        for k = 1:rho
            CDFs(k) <= CDFs(k+1);  % monotonic non-decreasing
        end
cvx_end

end