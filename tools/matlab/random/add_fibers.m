function [C, S, O] = add_fibers(r, n, m)

  rr = 1.3;

  x1 = -2;
  y1 = 2;
  x2 = 2;
  y2 = -2;

  lx = x2 - x1;
  ly = y1 - y2;

  F = [x1 + lx / 2, y2 + ly / 2];

  max_t = 1000;

  for i = 2 : n
    for j = 1 : max_t
      x = rand() * lx + x1;
      y = rand() * ly + y2;
      [~, b] = knnsearch(F, [x y]);
      if b > 2 * r * rr
        F(i,:) = [x y];
        break
      end
    end
  end

  C = [];
  S = [];
  O = [];
  cc = knnsearch(F, [0 0], 'K', m);
  ss = knnsearch(F, [0 0], 'K', 5 * m);
  for i = 1 : n
    if ismember(i, cc)
      C = [C; F(i,:)];
    elseif ismember(i, ss)
      S = [S; F(i,:)];
    else
      O = [O; F(i,:)];
    end
  end
end
