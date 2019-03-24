function [inner, outer] = edge_nodes(C, S, r)

  n = 500;
  d = 2 * pi / n;

  inner = zeros(size(C,1),2);
  outer = zeros(size(S,1),2);

  for i = 1 : size(C,1)
    for j = 1 : n
      inner((i-1)*n + j,:) = [r * cos(d * j), r * sin(d * j)] + C(i, :);
    end
  end

  for i = 1 : size(S,1)
    for j = 1 : n
      outer((i-1)*n + j,:) = [r * cos(d * j), r * sin(d * j)] + S(i, :);
    end
  end
end
