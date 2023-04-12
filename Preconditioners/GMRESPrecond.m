N = input('Enter dimensions N');  %N = 50;                    %dimensions
alpha = input('Enter alpha'); %alpha = 0.01;              %alpha
dx = input('Enter dx'); %dx = 0.1 ;                 %dx
A = A_Matrix(N,alpha,dx);    %A Matrix
for i = 1:N
    b(i) = sin(i*dx);        %setting b as sin(i*dx)
end
b = b';                      %or
%b = zeros(N,1);
%b(randi([1,N],1)) = 1;        %setting b as random step function (basis vector)

% Set up preconditioner M
M = diag(diag(A));

x{1} = zeros(N,1);            %x_0 = 0 , solution vectors
r{1} = b;                     %r_0 = b - Ax_0 , residual vectors
p{1} = r{1};                  %p_0 = r_0 , direction vectors

iter = 30        %set iterations for GMRES

for n = 1:iter
    alpha_val(n+1) = (p{n}' * A' * r{n}) / (p{n}' * A' * A * p{n});
    x{n+1} = x{n} + alpha_val(n+1)*p{n};          %updating x using alpha
    r{n+1} = r{n} - alpha_val(n+1)*A*p{n};        %updating r using alpha
    
    % Apply preconditioner to the residual
    z = M \ r{n+1};
    
    u = A*p{n};       %2nd krylov basis vector
    v = A*u;          %3rd krylov basis vector
    for j = 1:n
        beta = -(p{j}' * A' * v)/(p{j}' * A' * A *p{j});
        u = u + beta*p{j};        %updating u using beta
        v = v + beta*A*p{j};      %updating v using beta
    end
    
    % Apply preconditioner to the direction vector
    p{n+1} = M \ u;
end

x_star = x{end}       %our numerical solution
x_pre = p{end}