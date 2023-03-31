function [A] = A_Matrix(n,alpha,dx)     %n = dimension of A
                                        %alpha = constant
A = zeros(n,n);                         % dx = gridspace resolution
for i = 2:n-1
    A(i,i-1) = 1/(dx*dx);            %setting values for A excluding
    A(i,i) = -alpha - 2/(dx*dx);     %A_1i and A_ni
    A(i,i+1) = 1/(dx*dx);
end
A(1,1) = -alpha -2/(dx*dx);
A(1,2) = 1/(dx*dx);                %setting values for the 1st row and
A(1,n) = 1/(dx*dx);                %nth row
A(n,1) = 1/(dx*dx);
A(n,n-1) = 1/(dx*dx);
A(n,n) = -alpha - 2/(dx*dx);
end