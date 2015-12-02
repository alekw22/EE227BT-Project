clear all
load('short_time.mat');


X_new_1 = X_res_1(1:2:end,:);
X_new_2 = X_res_1(2:2:end,:);

batch = size(X_res_1,2);

violations  = zeros(T,1);

for j =1:T 
    
%calculate percent violation
x1_set = find(X_new_1(j+1,:) < x_tar(1,:));
x1_set = X_new_1(j+1,x1_set);

x2_set = find(X_new_2(j+1,:) < x_tar(2,:));
x2_set = X_new_2(j+1,x2_set);

intrsct_x = intersect(x1_set,x2_set);
union_x = union(x1_set,x2_set);

violation_temp = size(union_x,1)- size(intrsct_x,1);
violation_temp = (violation_temp*100)/(size(X_new_1(j,:),2));

violations(j,1) = violation_temp;
end

cost = sum(diag(X_new_1(2:end,:)'*X_new_1(2:end,:))) +sum(diag(U'*U));
cost = cost/batch;
