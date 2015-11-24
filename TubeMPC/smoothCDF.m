function CDFs = smoothCDF(CDF)
% smooth a discretized approximation of a Cumulative Discretization
% Function subject to constraints that make it a CDF

rho = length(CDF) - 1;

cvx_begin quiet
    variable CDFs(rho+1)
    minimize(norm(CDFs - CDF))     % curve fit
    subject to
    0 <= CDFs <= 1                 % must be a probability
        for k = 1:rho
            CDFs(k) <= CDFs(k+1);  % monotonic non-decreasing
        end
cvx_end

end