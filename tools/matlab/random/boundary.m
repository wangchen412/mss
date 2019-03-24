function p = boundary(theta, n, x)

  u = linspace(-x, x, 400);
  v = linspace(-x, x, 400);

  z = zeros(length(u), length(v));
  for i = 1:length(u)
    for j = 1:length(v)
      z(i,j) = map_feature(u(i), v(j), n)*theta;
    end
  end

  z = z';

  p = contour(u, v, z, [0, 0], 'LineWidth', 1);
  p = p(:, 2:size(p,2)-1);
  p = p';

end

