clc
clear all
z1 = read_complex("z1");
z2 = read_complex("z2");
[z2m, z1m] = precond(z2, z1);
check_ev(eig(z2m' * z2m, z2m' * z1m))