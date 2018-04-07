clc
clear all
z1 = read_complex("z1");
z2 = read_complex("z2");
check_ev(eig(z1' * z2, z1' * z1))