clc; close; clear all;

h = 8.1; % m
r = 2; % m
m = 9600; % kg

Ic_spaceship = [1/12*m*(3*r^2 + h^2) 0 0;
               0 1/12*m*(3*r^2 + h^2) 0;
               0 0 1/2*m*r^2]