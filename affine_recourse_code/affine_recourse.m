A = [0.7 -0.1*2;-0.1*3 0.9];
B = eye(2);
%five time steps
sizeA = size(A);
sizeAB = size(A*B);
I = eye(size(A));
G = [repmat(zeros(sizeA),1,5);
     I repmat(zeros(sizeA),1,4);
     A I repmat(zeros(sizeA),1,3);
     A^2 A I repmat(zeros(sizeA),1,2);
     A^3 A^2 A  I repmat(zeros(sizeA),1,1);
     A^4 A^3 A^2 A I];
 
 H = [repmat(zeros(sizeAB),1,5);
      B repmat(zeros(sizeAB),1,4);
      A*B B repmat(zeros(sizeAB),1,3);
      A^2*B A*B B repmat(zeros(sizeAB),1,2);
      A^3*B A^2*B A*B  B repmat(zeros(sizeAB),1,1);
      A^4*B A^3*B A^2*B A*B *B A^5*B];
  
 Sigma_little = 0.1*eye(2);
 
 x_d = ones(12,1);
 x_o = ones(2,1);
 x_zero = [x_o; A*x_o; A^2*x_o; A^3*x_o;A^4*x_o;A^5*x_o];
 
 w_high = ones(10,1);
 w_low = -w_high;
 
 Sigma = blkdiag(Sigma_little,Sigma_little,Sigma_little,Sigma_little,Sigma_little);
 Sigma = Sigma.^(1/2);
 n = 12;
 
 cvx_begin
 variable x(n);
 variable x_bar(n);
 variable u(n-2);
 variable u_bar(n-2);
 variable Q(n-2,n) lower_triangular;
 variable r(n-2);
 
 minimize((sum_square(x)) + (sum_square(u)))
 
 subject to
 
 %[G + H*Q*G; Q*G]*3*Sigma*w_low <= [x;u] - [x_bar;u_bar] <= [G + H*Q*G; Q*G]*3*Sigma*w_high;
 [x;u] -[x_bar;u_bar] == 0;
 u_bar == Q*x_zero + r;
 x_bar == (eye(size(H*Q)) + H*Q)*x_zero + H*r;
 x >= ones(12,1);
 -5*ones(10,1)<= u <= 5*ones(10,1);

 
 cvx_end
 
 K = ((eye(size(Q*H))+Q*H)^(-1))\Q ;
 u_o = (((eye(size(Q*H))+Q*H)^(-1)))\r;
 
 %simulate dynamics
 w = mvnrnd(zeros(10,1),Sigma.^2,1);
 
 x = (eye(size(H*K))- (H*K))^(-1)*(G*w' + H*u_o +x_zero);
 %%
 x_one = x([1:2:end],:);
 x_two = x([2:2:end],:);
 
 plot(x_one,x_two);
 