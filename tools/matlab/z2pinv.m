clc
clear all
z1 = read_complex("z1");
z2 = read_complex("z2");
z1g = read_complex("z1gp");
z2g = read_complex("z2gp");

[z1m, z2m] = precond(z1, z2);
check_ev(eig(pinv(z1m)*z2m))
check_ev(eig(pinv(z1g)*z2g))