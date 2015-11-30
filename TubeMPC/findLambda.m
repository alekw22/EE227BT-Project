clear; close all;

f = @(lambda) calcTubeParam(lambda);
[lambda2 y] = Bisection(0.2,0.99,5,f)

%lambda = 0.2741;
V2 = calcTubeParamV(lambda1)