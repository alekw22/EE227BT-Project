function [xmin,fmin]=Bisection(a,b,n,f)
%Performs the bisection algorithm


% First plot the function
fplot(f,[a,b],5);

hold on

% To display results

disp('Iterations        x1                 fx1          Df');
disp('-------------------------------------------------------');

% Start Iterations

k=1;

% If the number of iterations are less than the max. iterations
while k<=n
    x1=a;
    x2=b;
    x3=(a+b)/2;
    
% To find derivative we need to find 'deltax' 


if abs(x3)>.01
    deltax=.01*abs(x2);
else
    deltax=.001;
end

% Find the value of Derivatives at 'x1','x2' and 'x3'

    Df1=(f(x1+deltax)-f(x1-deltax))/(2*deltax);
    Df2=(f(x2+deltax)-f(x2-deltax))/(2*deltax);
    Df3=(f(x3+deltax)-f(x3-deltax))/(2*deltax);
    
% Region Elimination 
   
    if Df1*Df2>0
        disp('No minima in (a,b)')
    elseif (Df1*Df3)>0
        a=x3;
    elseif (Df2*Df3)>0
        b=x3;
    end

% Find 'x' at which 'f' is minimum

    xmin=x3;
    fmin=f(x3);
    
    
% Plot points

    plot(x3,f(x3),'or','MarkerSize',10);
    hold on
   
    % This will display the iterations on screen
    fprintf('%4d  %20e %14e %14e\n', k-1,x3,f(x3),Df3)
    k=k+1;
end