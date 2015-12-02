clear
T = 6;
epoch = 1;
eta = 3;
batch = 200;
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

X_res_1 = zeros(n*(epoch*T + 1), batch);
X_res_2 = zeros(n*(epoch*T + 1), batch);

U = zeros(n*epoch*T,batch);

for j = 1:batch
    x_01 = [1,1]';
    x_02 = x_01;
    X_1 = x_01;
    X_2 = x_02;
    x_1 = [];
    x_2 = [];
    U_res = [];
    status = '';
    
    for i = 1:epoch
        if i == 1
            
        else
            x_01 = x_1(end-n+1:end);
            x_02 = x_2(end-n+1:end);
            
        end
        
        while ~strcmp(status,'Solved')
            [F,u_0,status] = affine_recourse4(A,B,x_01,x_tar,n,T,mu,Sigma,eta*(0.1^0.5));
        end
        [x_1,x_2,u_res] = propagate_dynamics(A,B,T,u,x_01,x_02,F,u_0,mu,Sigma);
        X_1 = [X_1;x_1(n+1:end,:)];
        X_2 = [X_2;x_2(n+1:end,:)];
        U_res = [U_res;u_res];
        status = '';
    end
    
    X_res_1(:,j) = X_1;
    X_res_2(:,j) = X_2;
    U(:,j) = U_res;
    
    
end
% %%
% %plot result
% X_new_1 = X_1(1:2:end,:);
% X_new_2 = X_1(2:2:end,:);
%
% X_nfb_1 = X_2(1:2:end,:);
% X_nfb_2 = X_2(2:2:end,:);
%
% %plot(X_new_1',X_new_2');
%
% for j = 1:size(X_new_1,2)
%
%     plot(X_new_1(:,j), X_new_2(:,j));
%     hold on;
%
% end
% hold off
%
% figure
% for j = 1:size(X_nfb_1,2)
%
%     plot(X_nfb_1(:,j), X_nfb_2(:,j));
%     hold on;
%
% end
% title('no feedback');
%
% %calculate percent violation
% x1_set = find(X_new_1 < x_tar(1,:));
% x1_set = X_new_1(x1_set,:);
%
% x2_set = find(X_new_2 < x_tar(2,:));
% x2_set = X_new_2(x2_set,:);
%
% intrsct_x = intersect(x1_set,x2_set);
% union_x = union(x1_set,x2_set);
%
% violation = size(union_x,1)- size(intrsct_x,1);
% violation = (violation*100)/(size(X_new_1,1))
%
% hold off
% figure
% hold on
%
% h =  animatedline;
%
%
%
% for j = 1:size(X_new_1,1)
%     addpoints(h,X_new_1(j,:),X_new_2(j,:));
%     pause(0.2);
%     drawnow;
% end
% x1 = get(gca,'xlim');
% plot(x1,[1,1]);
%
% y1 = get(gca,'ylim');
% plot([1,1],y1);
% hold off;
% %%
% figure;
% hold on
%
% b = animatedline;
% for j = 1:size(X_nfb_1,1)
%     addpoints(b,X_nfb_1(j,:),X_nfb_2(j,:));
%     pause(0.2);
%     drawnow;
% end
%
% x1 = get(gca,'xlim');
% plot(x1,[1,1]);
%
% y1 = get(gca,'ylim');
% plot([1,1],y1);
%
