clc
clear all
z1 = read_complex("z1");
z2 = read_complex("z2");
[z1m, z2m] = precond(z1, z2);
check_ev(eig(z1m' * z2m, z1m' * z1m))