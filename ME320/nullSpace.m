clear all; clc;
syms a b c
A = [0,1,1;
    1,-1,-2;
    0,-1,-2];
w = [0;-3;1];
x = [a,b,c];
NullSpace = null(A,'r');
linsolve(A,w)
