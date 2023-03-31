function [] = GMRES(A,b)
%GMRES Summary of this function goes here
%   Detailed explanation goes here

N = size(A,1);%calcs size

alpha = zeros(1,N+1);%sets alpha, x, r, and p
x = cell(1,N+1);
r = cell(1,N+1);
p = cell(1,N+1);


x{1} = zeros(N,1);%initialises x
r{1} = b';
p{1} = r{1};

for n = 1:N%loops N times
    alpha(n+1) = (p{n}' * A * r{n}) / (p{n}' * A * p{n});%finds new alpha, x and r
    x{n+1} = x{n} - alpha(n+1)*p{n};
    r{n+1} = r{n} - alpha(n+1)*A*p{n};
    u = A*p{n};
    v = A*u;
    for j = 1:n%loops again n times
        beta = -(p{j}'*A*v)/(p{j}'*A'*A*p{j});%sets beta, u, and v
        u = u + beta*p{j};
        v = v + beta*A*p{j};
    end
    p{n+1} = u;%sets p
end

end

