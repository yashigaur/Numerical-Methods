clear all;
clc;

g = 9.81; %acceleration due to gravity (m/s^2)
L = 30; %length (m)
m = 68.1; %mass of person (kg)
cd = 0.25; %drag coefficient (kg/m)
k = 40; %spring constant (N/m)
p = 8; %cord dampening constant (kg/s)

n=500 ; %total time steps

dt = 0.1; %timestep 

%-------------------------------------------------------------------------%
%Euler's method
time1 = zeros(n,1) ; %time vector
x1= zeros(n,1); %displacement vector
v1= zeros(n,1); %velocity vector

%initial conditions
x1(1)=0;
v1(1)=0;
time1(1)=0;

for i=1:n
    time1(i+1)=time1(i)+dt;
    t=0;
    if(x1(i)>L)
        t = (k/m)*(x1(i)-L)+(p/m)*v1(i);
    end
    v1(i+1) = v1(i) + (g - sign(v1(i))*(cd/m)*(v1(i)^2) - t)*dt;
    x1(i+1) = x1(i) + v1(i)*dt;
end

figure;
plot(time1,v1,'DisplayName','velocity');
hold on;
plot(time1,-x1,'DisplayName','position');
legend
title("Euler's method");
xlabel("time");

%-------------------------------------------------------------------------%
%Heun's method
time2 = zeros(n,1) ; %time vector
x2= zeros(n,1); %displacement vector
v2= zeros(n,1); %velocity vector

%initial conditions
x2(1)=0;
v2(1)=0;
time2(1)=0;

for i=1:n
     time2(i+1)=time2(i)+dt;
     t=0;
     x2(i+1) = x2(i) + dt;
     if(x2(i)>L)
         t = (k/m)*(x2(i)-L)+(p/m)*v2(i);
     end
     r =(g - sign(v2(i))*(cd/m)*(v2(i)^2) - t);   %slope 1
     v2(i+1) = v2(i) + (g - sign(v2(i))*(cd/m)*(v2(i)^2) - t)*dt; %predictor solution
     t0=0;
     if(x2(i)>L)
         t0 = (k/m)*(x2(i+1)-L)+(p/m)*v2(i+1);
     end
     p = (g - sign(v2(i+1))*(cd/m)*(v2(i+1)^2) - t0);  %slope 2
     v2(i+1) = v2(i) + 0.5*(r+p)*dt;  %corrector solution
     x2(i+1) = x2(i) + v2(i)*dt;

end
figure;
plot(time2,v2,'DisplayName','velocity');
hold on;
plot(time2,-x2,'DisplayName','position');
legend
title("Heun's method");
xlabel("time");

%-------------------------------------------------------------------------%
%Midpoint method
time3 = zeros(n,1) ; %time vector
x3= zeros(n,1); %displacement vector
v3= zeros(n,1); %velocity vector

%initial conditions
x3(1)=0;
v3(1)=0;
time3(1)=0;

for i=1:n
    time3(i+1)=time3(i)+dt;
    x3(i+1) = x3(i) + v3(i)*dt;
    a = x3(i) + dt/2;
    b = v3(i) + dt/2;
    q=0;
    if(a>L)
        q = (k/m)*(a-L)+(p/m)*b;
    end
    v3(i+1) = v3(i) + (g - sign(b)*(cd/m)*(b^2) - q)*dt;  %corrector solution

end
figure;
plot(time3,v3,'DisplayName','velocity');
hold on;
plot(time3,-x3,'DisplayName','position');
legend
title("midpoint method");
xlabel("time");

%-------------------------------------------------------------------------%
%Range kutta method
time4 = zeros(n,1) ; %time vector
x4= zeros(n,1); %displacement vector
v4= zeros(n,1); %velocity vector

%initial conditions
x4(1)=0;
v4(1)=0;
time4(1)=0;
%for this question i created a special function called slope, that takes
%two varibales x and v and returns me the value of dv/dt for those values
%function can be found in slope.m
for i=1:n
    time4(i+1)=time4(i)+dt;
    x4(i+1) = x4(i) + v4(i)*dt;
    k1 = slope(x4(i),v4(i));
    k2 = slope(x4(i) + dt/2,v4(i)+(dt/2)*k1);
    k3 = slope(x4(i) + dt/2,v4(i)+(dt/2)*k2);
    k4 = slope(x4(i) + dt,v4(i)+(dt)*k3);
    v4(i+1) = v4(i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end
figure;
plot(time4,v4,'DisplayName','velocity');
hold on;
plot(time4,-x4,'DisplayName','position');
legend
title("4th order Range Kutta method");
xlabel("time");

%plotting all plots in one graph
figure;
subplot(2,1,1);
plot(time1,v1);
hold on;
plot(time2,v2);
plot(time3,v3);
plot(time4,v4);
hold off;
legend({'Euler','Heun','Midpoint','Range-Kutta'})
title('Velocity using 4 methods');

subplot(2,1,2);
plot(time1,-x1);
hold on;
plot(time2,-x2);
plot(time3,-x3);
plot(time4,-x4);
hold off;
legend({'Euler','Heun','Midpoint','Range-Kutta'})
title('Position using 4 methods');
