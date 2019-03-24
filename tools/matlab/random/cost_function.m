function [J, grad] = cost_function(theta, X, y)
  h = 1 ./ (1 + exp(-X * theta));
  J = -(y' * log(h) + (1 - y)' * log(1 - h)) / length(y);
  grad = X' * (h - y) / length(y);
end
