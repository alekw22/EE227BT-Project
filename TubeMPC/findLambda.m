clear; close all;

f = @(lambda) calcTubeParam(lambda);
fplot(f,[0.2,0.99],5);

[lambda1 y] = Bisection(0.9,0.99,5,f)

V1 = calcTubeParamV(lambda1)