function out = map_feature(X1, X2, n)
  out = ones(size(X1(:,1)));
  for i = 1 : n
    for j = 0:i
      out(:, end+1) = (X1.^(i-j)).*(X2.^j);
    end
  end
end
