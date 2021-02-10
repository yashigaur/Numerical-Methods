clear all; % clear everything
clc; % clear command window
% define the parameters
To=21; % initial tempertaue
k=0.017; % proportionality constant

% dependent and independent variables
n=11; % number of time steps
T=zeros(n,1); % array for temperature computed numerically
Texact=zeros(n,1); % array of exact temperature
t=zeros(n,1);

% Initial condictions
T(1)=68;
Texact(1)=68;
t(1)=0;

dt=1; % seconds
%t(2)=t(1)+dt;

% time marching and solution computation
for i=2:n
    t(i)=t(i-1)+dt;
    Texact(i)=To+ 47*exp(-k*t(i)) ;
    T(i)=T(i-1)- k*(T(i-1)-To)*dt;
    end

% plotting
plot(t,Texact);
hold on;
plot(t,T);

Temperature.m
Displaying Temperature.m.