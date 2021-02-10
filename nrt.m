clear all; 
clc; 
% define the parameters
m=90; % mass in kg
h=0.45; % height in m
g=9.81; % accelaration due to gravity in m/s2
k1=50000; % force constant 1
k2=40; %force constant 2


xg=0.1; %guess value
fxg=((2*k2*(xg.^2.5))/5)+((k1*(xg.^2))/2)-(m*g*(xg))-m*g*h;

while ((abs(fxg))>0.001)
    fxg=((2*k2*(xg.^2.5))/5)+((k1*(xg.^2))/2)-(m*g*(xg))-m*g*h;
    fdashxg=(k2*(xg.^1.5))+(k1*xg)-(m*g);
    xg=xg-(fxg/fdashxg);
end



        