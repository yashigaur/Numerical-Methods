clear all; 
clc; 
% define the parameters
m=90; % mass in kg
h=0.45; % height in m
g=9.81; % accelaration due to gravity in m/s2
k1=50000; % force constant 1
k2=40; %force constant 2

%decide the bracket
xl=0.1;
xu=0.3;
xr=(xl+xu)/2;
fxr=((2*k2*(xr.^2.5))/5)+((k1*(xr.^2))/2)-(m*g*(xr))-m*g*h;

%bisection algorithm
while   (abs(fxr)>0.01)
    fxr=((2*k2*(xr.^2.5))/5)+((k1*(xr.^2))/2)-(m*g*(xr))-m*g*h;
    fxl=((2*k2*(xl.^2.5))/5)+((k1*(xl.^2))/2)-(m*g*(xl))-m*g*h;
    if ((fxl*fxr)==0)
        break;
    end
    if ((fxl*fxr)<0)
        xu=xr;
        xr=(xl+xr)/2;
    else
        xl=xr;
        xr=(xl+xu)/2;
    end
end
xr



        