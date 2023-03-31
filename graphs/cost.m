function [J] = cost(y_hat,y)

J = [];
n = length(y_hat);
    J = (1/n) * sum((y_hat - y).^2);
end
