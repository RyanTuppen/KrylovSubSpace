clear all;

N = input('Enter dimensions N');  %N = 50;                    %dimensions
alpha = input('Enter alpha'); %alpha = 0.01;              %alpha
dx = input('Enter dx'); %dx = 0.1 ;                 %dx
A = A_Matrix(N,alpha,dx);    %A Matrix
%for i = 1:N
%    b(i) = sin(i*dx);        %setting b as sin(i*dx)
%end
%b = b';                      %or
b = zeros(N,1);
b(randi([1,N],1)) = 1;        %setting b as random step function (basis vector)
x{1} = zeros(N,1);            %x_0 = 0 , solution vectors
r{1} = b;                     %r_0 = b - Ax_0 , residual vectors
p{1} = r{1};                  %p_0 = r_0 , direction vectors

iter = 30        %set iterations for GMRES

for n = 1:iter
    alpha_val(n+1) = (p{n}' * A' * r{n}) / (p{n}' * A' * A * p{n});
    x{n+1} = x{n} + alpha_val(n+1)*p{n};          %updating x using alpha
    r{n+1} = r{n} - alpha_val(n+1)*A*p{n};        %updating r using alpha
    u = A*p{n};       %2nd krylov basis vector
    v = A*u;          %3rd krylov basis vector
    for j = 1:n
        beta = -(p{j}' * A' * v)/(p{j}' * A' * A *p{j});
        u = u + beta*p{j};        %updating u using beta
        v = v + beta*A*p{j};      %updating v using beta
    end
    p{n+1} = u;
end

x_star = x{end}       %our numerical solution

%cost function
J = [];
for i = 1:iter + 1%iterates until converged
    J(i) = cost(A*x{i},b);       %calculating J / MSE values
end

iterations= [0:iter];%initialises iterations (x-axis)
figure (2);                      %plot
plot(iterations,J);
xlabel('Number of Iterations');%labels axis
ylabel('J');
title('Mean Squared Error against iteration for GMRES method - Varying alpha')%gices title
hold on;%allows multiple plots
legend('alpha = 0.5', 'alpha = 1', 'alpha = 2',  'alpha = 3', 'alpha = 5');%adds legend
exportgraphics(gca,"GMRES_varying_resolution.jpg");%saves image