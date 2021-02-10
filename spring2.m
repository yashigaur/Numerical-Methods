clear all; % clear everything
clc; % clear command window
% define the parameters
m=0.09; % mass in kg
h=0.45; % height in m
g=9.81; % accelaration due to gravity in m/s2
k1=50; % force constant 1
k2=0.04; %force constant 2

d=linspace(0,1,10);
f=((2*k2*(d.^2.5))/5)+((k1*(d.^2))/2)-(m*g*(d))-m*g*h;

plot(d,f);

spring2.m
Displaying spring2.m.

