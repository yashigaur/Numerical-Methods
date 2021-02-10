clear all;
clc;


n=500; %timestep
%initialising matrices
time = zeros(n+1,1);
c1 = zeros(n+1,1); 
c2 = zeros(n+1,1);
c3 = zeros(n+1,1);

%defining initial values
time(1)=0;
c1(1)=1;
c2(1)=1;
c3(1)=0;

dt=0.1;
%using implicit Euler's method
for i=1:n+1
    time(i+1)= time(i)+dt;
    c1(i+1) = c1(i)/(1 + 0.013*dt + 1000*dt*c3(i));  %defining C1 
    c2(i+1) = c2(i)/(1 + 2500*c3(i)*dt); %defining C2
    c3(i+1) = (c3(i) - 0.013*c1(i)*dt)/(1 + 1000*c1(i)*dt + 2500*c2(i)*dt); %defining C3
end

%plotting
plot(time,c1);
hold on;
plot(time,c2);
plot(time,c3,"");
xlabel("Time (s)");
ylabel("Concentration (c)");
ylim([-0.2,2]);
legend("C1","C2","C3");