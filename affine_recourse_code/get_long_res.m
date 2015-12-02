clear all
load('long_time.mat');

X_new_1 = X_res_1(1:2:end,:);
X_new_2 = X_res_1(2:2:end,:);

%calculate percent violation
x1_set = find(X_new_1(2:end,:) < x_tar(1,:));
x1_set = X_new_1(x1_set,:);

x2_set = find(X_new_2(2:end,:) < x_tar(2,:));
x2_set = X_new_2(x2_set,:);

intrsct_x = intersect(x1_set,x2_set);
union_x = union(x1_set,x2_set);

violation = size(union_x,1)- size(intrsct_x,1);
violation = (violation*100)/(size(X_new_1,1)-1);

cost = norm(U(1:end,1))^2 + norm(X_res_1(1:end-2,1))^2;

figure
hold on;
xlabel('x1');
ylabel('x2');
title('x1 vs x2 affine recourse');


for j = 1:size(X_new_1,2)

    plot(X_new_1(:,j), X_new_2(:,j),'.', 'MarkerSize', 7);
    

end


xlim([0,4])

ylim([0,4])

plot([1,1],[0 4], 'LineWidth', 1, 'Color', 'k')

plot([0 4],[1,1], 'LineWidth', 1, 'Color', 'k')

hold off;
