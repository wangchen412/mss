clc
clear

rr = 1.3;
r = 0.06;
x1 = -4;
y1 = 4;
x2 = 4;
y2 = -4;
lx = x2 - x1;
ly = y1 - y2;

F = [x1 + lx / 2, y2 + ly / 2];
max_t = 1000;

n = 1600;

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

min_dist = 10;

for i = 1 : n
  [~, b] = knnsearch(F, F(i,:), 'K', 2);
  if b(2) < min_dist
    min_dist = b(2);
  end
end

F = [F, zeros(size(F, 1), 1)];