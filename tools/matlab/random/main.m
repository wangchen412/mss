clear all
close all
clc

n = 6;

[C, S, O] = add_fibers(0.06, 400, 40);
[inner, outer] = edge_nodes(C, S, 0.06);
X = [inner; outer];
y = [ones(size(inner,1),1);zeros(size(outer,1),1)];
X = map_feature(X(:,1), X(:,2), n);
options = optimset('Algorithm','trust-region','GradObj', 'on', 'MaxIter', 100000);
[theta, J, exit_flag] = fminunc(@(t)(cost_function(t, X, y)), zeros(size(X, 2), 1), options);

plot_fibers(C, S, O)
hold on
boundary(theta, n, 1);

figure();
hold on
scatter(inner(:,1), inner(:,2), 1, "r");
scatter(outer(:,1), outer(:,2), 1, "b");
p = boundary(theta, n, 1);
F = [C; S; O];
