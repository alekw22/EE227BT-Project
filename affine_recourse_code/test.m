clear
T = 10;
epoch = 10;
n = 2;
u = zeros(n*T,1);
A = [0.7 -0.1*2;-0.1*3 0.9];
%A = [0,-1;1,0];
B = eye(n);
x_01 = [1,1]';
x_02 = x_01;
x_tar = [1,1]';
mu = zeros(2,1);
Sigma = 0.1*eye(2,2);
X_1 = x_01;
X_2 = x_02;
x_1 = [];
x_2 = [];
status = '';
for i = 1:epoch
    if i == 1
        
    else
        x_01 = x_1(end-n+1:end);
        x_02 = x_2(end-n+1:end);
        
    end
    
    while ~strcmp(status,'Solved')
        [F,u_0,status] = affine_recourse4(A,B,x_01,x_tar,n,T,mu,3*Sigma,3*(0.1^0.5));
    end
    [x_1,x_2] = propagate_dynamics(A,B,T,u,x_01,x_02,F,u_0,mu,Sigma);
    X_1 = [X_1;x_1(n+1:end,:)];
    X_2 = [X_2;x_2(n+1:end,:)];
    status = '';
end

%plot result
X_new_1 = X_1(1:2:end,:);
X_new_2 = X_1(2:2:end,:);

X_nfb_1 = X_2(1:2:end,:);
X_nfb_2 = X_2(2:2:end,:);

%plot(X_new_1',X_new_2');

for j = 1:size(X_new_1,2)
    
    plot(X_new_1(:,j), X_new_2(:,j));
    hold on;
    
end
hold off

figure
for j = 1:size(X_nfb_1,2)
    
    plot(X_nfb_1(:,j), X_nfb_2(:,j));
    hold on;
    
end
title('no feedback');

%calculate percent violation
x1_set = find(X_new_1 < x_tar(1,:));
x1_set = X_new_1(x1_set,:);

x2_set = find(X_new_2 < x_tar(2,:));
x2_set = X_new_2(x2_set,:);

intrsct_x = intersect(x1_set,x2_set);
union_x = union(x1_set,x2_set);

violation = size(union_x,1)- size(intrsct_x,1);
violation = (violation*100)/(size(X_new_1,1))

