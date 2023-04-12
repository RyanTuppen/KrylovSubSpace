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

M = diag(diag(A)); %preconditioner, diagonal of A
Minv = inv(M); %inverse of preconditioner

x{1} = zeros(n,1);          %x_1 set as zero vector where x = solution iterations
r{1} = b - A*x{1};          %r_1 set as b, where r = residuals
z{1} = Minv * r{1};         %apply preconditioner to initial residual
p{1} = z{1};                %p_1 set as z_1, where p = directions

alpha(1) = 0;
beta(1) = 0;

N = 30;                  %30 iterations
for i = 1:N
    alpha(i+1) = (r{i}' * z{i}) / (p{i}' * A * p{i});
    x{i+1} = x{i} + alpha(i+1) * p{i};                     %updating our solution and
    r{i+1} = r{i} - alpha(i+1) * A * p{i};                 %residual using alpha
    z{i+1} = Minv * r{i+1};                                %apply preconditioner to residual
    beta(i+1) = (z{i+1}' * r{i+1}) / (z{i}' * r{i});       %compute beta using preconditioned vectors
    p{i+1} = z{i+1} + beta(i+1) * p{i};                %updating our
end                                                    %direction using beta

x{end}