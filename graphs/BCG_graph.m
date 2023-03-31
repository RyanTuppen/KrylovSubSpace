%Bi-conjugate method plotting graph
clear all
%n = input('Enter dimensions N');  %n = 10;                    %dimensions
%alpha = input('Enter alpha'); %alpha = 2;              %alpha
%dx = input('Enter dx'); %dx = 0.1 ;                 %dx

n = 50;
alpha = 2;   
dx = 0.3; 

A = A_Matrix(n,alpha,dx);       %A matrix          

for i = 1:n
   b(i) = sin(i*dx);        %known vector b set as sin(x_i) or random
end                         %basis vector
b = b';

%b = zeros(n,1);
%b(randi([1,n],1)) = 1;

res_tol  = 1e-9;%sets resolution tolerance
max = 1000; %sets maximum number of iterations
N = size(A, 1);%Finds size of A
x = zeros(N, 1);%sets vector x

r = b - A * x;%calculates r
p = r;%sets values
rs = r;
ps = p;
r2 = r' * rs;%sets r2 so can be interchanged later

J=[];
i=1;
	
residual = norm(r, 2);%calculates residual norm
res_stop = residual * res_tol;%stop condition
count = 1;%initialises count
res_norm(count) = residual;%sets reolution norm
	
converged = 0;%initialises converged variable
while ((count < max) && (residual > res_stop))%loops until stop conditions met
	alpha = r2 / ((A*p)' * ps);%Calculates alpha
	x = x + alpha * p;%adds onto x
	r = r - alpha * (A*p);%Finds r
	rs = rs - alpha * (A' * ps);%Finds rs

	%sets new and old r2
	r2_old = r2;
	r2 = r' * rs;
	beta = r2 / r2_old;
		
    %finds p and ps
	p  = r  + beta * p;
	ps = rs + beta * ps;
		
	count = count + 1;%increments count
	residual = norm(r, 2);%calculates new residual
	res_norm(count) = residual;

    J(i) = cost(A*x, b);%finds cost function each iteration
    i = i+1;%increments i
end

if residual <= res_stop %sets converged to 1 if residual less than stop
    converged = 1; 
end

iter=[0:i-2];
figure(1);                %plot
plot(iter,J);
xlabel('Number of iterations');
ylabel('J');
title('Mean Squared Error against Iterations for Bi-CG Method - Varying resolution (dx)');
hold on;%allows multiple plots
legend('dx=0.01', 'dx = 0.05', 'dx = 0.1',  'dx = 0.2');%adds legend
exportgraphics(gca,"BCG_varying_alpha_sin.jpg");%saves image
