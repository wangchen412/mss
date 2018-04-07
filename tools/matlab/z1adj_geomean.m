clc
clear all
z1g = read_complex("z1gp");
z2g = read_complex("z2gp");
check_ev(eig(z1g' * z2g, z1g' * z1g))