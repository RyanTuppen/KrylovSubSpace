clear all
n = input('Enter dimensions N');  %n = 10;                    %dimensions
alpha = input('Enter alpha'); %alpha = 2;              %alpha
dx = input('Enter dx'); %dx = 0.1 ;                 %dx
A = A_Matrix(n,alpha,dx);       %A matrix          

for i = 1:n
   b(i) = sin(i*dx);        %known vector b set as sin(x_i) or random
end                         %basis vector
b = b';

%b = zeros(n,1);
%b(randi([1,n],1)) = 1;

x{1} = zeros(n,1);          %x_1 set as zero vector where x = solution iterations
r{1} = b - A*x{1};          %r_1 set as b, where r = residuals
p{1} = r{1};                %p_1 set as r_1, where p = directions

alpha(1) = 0;
beta(1) = 0;

N = 30;        %set iterations for CG
for i = 1:N
    alpha(i+1) = (r{i}' * r{i}) / (p{i}' * A * p{i});
    x{i+1} = x{i} + alpha(i+1) * p{i};                     %updating our solution and
    r{i+1} = r{i} - alpha(i+1) * A * p{i};                 %residual using alpha
    beta(i+1) = (r{i+1}' * r{i+1}) / (r{i}' * r{i});
    p{i+1} = r{i+1} + beta(i+1) * p{i};                %updating our
end                                                    %direction using beta

x_star = x{end}

%cost function

J = [];
for i = 1:N+1
    J(i) = cost(A*x{i},b);       %calculating J / MSE values
end

iter= [0:N];
figure(1);                %plot
plot(iter,J);
xlabel('Number of iterations');
ylabel('J');
title('Mean Squared Error against Iterations for CG Method - Varying resolution (dx)');
hold on;%allows multiple plots
legend('dx = 0.05', 'dx = 0.1',  'dx = 0.2', 'dx = 0.3', 'dx = 0.4');%adds legend
%exportgraphics(gca,"CG_varying_resolution.jpg");%saves image

