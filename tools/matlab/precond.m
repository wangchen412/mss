function [z1_p, z2_p] = precond(z1, z2)
    z1_p = zeros(size(z1));
    z2_p = zeros(size(z2));
    for i = 1 : size(z1, 1)
        p = mean(z1(i, :));
        z1_p(i, :) = z1(i, :) ./ p;
        z2_p(i, :) = z2(i, :) ./ p;
    end
end

