%Conjugate Gradient method

n = input('Enter dimensions N');  %n = 10;                    %dimensions
alpha = input('Enter alpha'); %alpha = 2;              %alpha
dx = input('Enter dx'); %dx = 0.1 ;                 %dx

%n = 10;                      %dimensions
%alpha = 2;                   %constant
%dx = 1;                      %gridspace resolution

A = A_Matrix(n,alpha,dx);       %A matrix
b = ones(n,1);              
for i = 1:n                  %known vector b
    b(i) = sin(i*dx);        %set as sin(x_i)
end
x{1} = zeros(n,1);          %x_1 set as zero vector where x = solution iterations
r{1} = b - A*x{1};          %r_1 set as b, where r = residuals
p{1} = r{1};                %p_1 set as r_1, where p = directions

alpha(1) = 0;
beta(1) = 0;

N = 30;                  %30 iterations
for i = 1:N
    alpha(i+1) = (r{i}' * r{i}) / (p{i}' * A * p{i});
    x{i+1} = x{i} + alpha(i+1) * p{i};                     %updating our solution and
    r{i+1} = r{i} - alpha(i+1) * A * p{i};                 %residual using alpha
    beta(i+1) = (r{i+1}' * r{i+1}) / (r{i}' * r{i});
    p{i+1} = r{i+1} + beta(i+1) * p{i};                %updating our
end                                                    %direction using beta

x{end}