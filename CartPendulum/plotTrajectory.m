%clear; close all

load resultEps1
%load resultEps2
load cart_pend_tube

figure(2)
hold on
plot(xCL(1,:))
plot(xCL(3,:))
plot(x_tube(1,:))
plot(x_tube(3,:))

legend('cart pos (SCMPC)', 'angle (SCMPC)', 'cart pos (Tube)', 'angle (Tube)')