clear; clf;

load('results3000.mat')
figure(1)
hold on
plot(x(1,:)',x(2,:)','.', 'MarkerSize', 15)
xlim([0,4])
ylim([0,4])
plot([1,1],[0 4], 'LineWidth', 1, 'Color', 'k')
plot([0 4],[1,1], 'LineWidth', 1, 'Color', 'k')


c = sum(sum(Cost))
viol_q1 = sum(sum(v1))/3000 *100

viol_q2 = sum(sum(v2))/3000
total = (sum(sum(v1)) + sum(sum(v2)))/6000