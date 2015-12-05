clear
T = 6;
epoch = 25;
eta = 1.281551565545;
batch = 200;
load('./CartPendulum/sysDyn.mat');
%A = [0.7 -0.1*2;-0.1*3 0.9];
%B = eye(2);
%A = [1.6,1.1;-0.7,1.2];
%B = [1;1];
%Bw = eye(2);
%Bw = B;
n = size(A,2);
nn = size(B,2);
nw = size(Bw,2);
%A = [0,-1;1,0];

u = zeros(nn*T,1);
%x_01 = [-5,60]';
x_01 = [1,0,0,0]';
%x_02 = x_01;
x_tar = [0,0,0,0]';
%x_tar = [0,0]';
mu = zeros(nw,1);
Sigma = 0.1*eye(nw,nw);

X_res_1 = zeros(n*(epoch*T + 1), batch);
X_res_2 = zeros(n*(epoch*T + 1), batch);

U = zeros(nn*epoch*T,batch);


for j = 1:batch
    %x_01 = [-5,60]';
    x_01 = [1,0,0,0]';
    x_02 = x_01;
    X_1 = x_01;
    X_2 = x_02;
    x_1 = [];
    x_2 = [];
    U_res = [];
    status = '';
    gamma =1;
    j
    for i = 1:epoch
        i;
       
        if i == 1
            
        else
            x_01 = x_1(end-n+1:end);
            x_02 = x_2(end-n+1:end);
            
        end
        
        while ~(strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved'))
            [F,u_0,status] = affine_recourse_cart(A,B,Bw,x_01,x_tar,n,T,mu,Sigma,eta*(1.281551565545)^0.5,gamma);
        end
        [x_1,x_2,u_res] = propagate_dynamics_cart(A,B,Bw,T,u,x_01,x_02,F,u_0,mu,Sigma);
        X_1 = [X_1;x_1(n+1:end,:)];
        X_2 = [X_2;x_2(n+1:end,:)];
        U_res = [U_res;u_res];
        status = '';
        gamma = gamma*1.000;
    end
    
    X_res_1(:,j) = X_1;
    X_res_2(:,j) = X_2;
    U(:,j) = U_res;
    
    
end
%%
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
