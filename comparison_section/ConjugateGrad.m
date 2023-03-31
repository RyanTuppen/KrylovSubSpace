function [] = ConjugateGrad(A,b)
%CONJUGATEGRAD Summary of this function goes here
%   Detailed explanation goes here


n = size(A,1);%calcs size of A

N = 30;                  %30 iterations

alpha = cell(1,N+1);%sets alpha
x = cell(1,N+1);
r = cell(1,N+1);
beta = cell(1,N+1);
p = cell(1,N+1);

x{1} = zeros(n,1);          %x_1 set as zero vector where x = solution iterations
r{1} = b - A*x{1};          %r_1 set as b, where r = residuals
p{1} = r{1};                %p_1 set as r_1, where p = directions
alpha{1} = 0;
beta{1} = 0;

for i = 1:N
    alpha{i+1} = (r{i}' * r{i}) / (p{i}' * A * p{i});
    x{i+1} = x{i} + alpha{i+1} * p{i};                     %updating our solution and
    r{i+1} = r{i} - alpha{i+1} * A * p{i};                 %residual using alpha
    beta{i+1} = (r{i+1}' * r{i+1}) / (r{i}' * r{i});
    p{i+1} = r{i+1} + beta{i+1} * p{i};                %updating our
end                                                    %direction using beta

x{end};

end

