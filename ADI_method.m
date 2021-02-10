%Soving complex differential equation using Explicit method and Alternating Direction Implicit (ADI) method.
clear all;
clc;

alpha = 0.835; %thermal diffusivitiy

L = 100; %length of square plate (in cm)
time = 500; %total time (sec)

nt = 4000;%timesteps
nx = 50;%steps in x
ny = 50;%steps in y

dt = time/nt;
dx = L/nx;
dy = L/ny;

%Stability criterion 
%dt = (dx^2 + dy^2)/(8*alpha); 
% solving this for dx=1 and dy=1 we get dt<=0.299
%therfore dt=0.25 is stable

gammax = (alpha*dt)/(dx^2);
gammay = (alpha*dt)/(dy^2);

%nt=round(time/dt); %no. of timesteps

%explicit method
Temp1 = zeros(nx+1,ny+1,nt+1);

%initial conditions
Temp1(:,1,:)= 0;
Temp1(:,ny+1,:)= 100;
Temp1(1,:,:)= 75;
Temp1(nx+1,:,:)= 50;

%for time loop
for n=1:nt+1
    for i=2:nx
        for j=2:ny
            Temp1(i,j,n+1) = Temp1(i,j,n) + gammax*(Temp1(i+1,j,n) - 2*(Temp1(i,j,n)) + Temp1(i-1,j,n)) + gammay*(Temp1(i,j+1,n) - 2*(Temp1(i,j,n)) + Temp1(i,j-1,n));
        end
    end
end

t = input("Enter the time at which you want to see temperature variation (1 - 500): ");
[M,N] = size(Temp1(:,:,4*t));
[x,y] = meshgrid(1:N,1:M); 

subplot(2,2,3);
contour(x,y,Temp1(:,:,4*t));
colorbar;
title("Contour plot")
subplot(2,2,2);
surf(x,y,Temp1(:,:,4*t));
colorbar;
title(["Temp variation at ",num2str(t),"secs"]);
subplot(2,2,4);
surf(x,y,Temp1(:,:,4*t));
view(0,15);
colorbar;
sgtitle("Using Explicit Method")

%--------------------------------------------------------------------------
%ADI method
Temp2 = zeros(nx+1,ny+1,nt+1);  %defining temperature matrix

%boundary conditions
Temp2(:,1,:)= 0;
Temp2(:,ny+1,:)= 100;
Temp2(1,:,:)= 75;
Temp2(nx+1,:,:)= 50;

for n=1:nt+1
    %defining T_half variable for column sweep
    T_half = zeros(nx+1,ny+1); 
    T_half(:,1)= 0;
    T_half(:,ny+1)= 100;
    T_half(1,:)= 75;
    T_half(nx+1,:)= 50;
    lambda2 = gammay/2;
    lambda1 = gammax/2;
    %coloumn sweep
    for i=2:nx
        %essential parameters for thomas algo
        main_diag = zeros(ny-1,1);
        upp_diag = zeros(ny-1,1);
        lower_diag = zeros(ny-1,1);
        main_diag(:)= (1+2*lambda2);
        upp_diag(1:ny-2) = -lambda2;
        lower_diag(2:ny-1) = -lambda2;
        
        %definig the Rhs matrix
        R = zeros(ny-1,1);
        R(1) = Temp2(i,2,n) + lambda1*(Temp2(i+1,2,n) - 2*(Temp2(i,2,n)) + Temp2(i-1,2,n)) + 0*lambda2;
        R(ny-1)= Temp2(i,ny,n) + lambda1*(Temp2(i+1,ny,n) - 2*(Temp2(i,ny,n)) + Temp2(i-1,ny,n)) + 100*lambda2;
        for j=2:ny-2
            R(j) = Temp2(i,j+1,n) + lambda1*(Temp2(i+1,j+1,n) - 2*(Temp2(i,j+1,n)) + Temp2(i-1,j+1,n));
        end
        %applying thomas algo
        X = thomas(upp_diag,lower_diag,main_diag,R,ny-1);
        
        %final value of t_half
        for k=2:ny
            T_half(i,k)=X(k-1);
        end
    end
    %row sweep
    %as the geometry is symmetric we can apply same methods for both rows
    %and columns
    
    for j=2:ny
        %intoducing parameters for thomas algo
        main_diag = zeros(nx-1,1);
        upp_diag = zeros(nx-1,1);
        lower_diag = zeros(nx-1,1);
        main_diag(:)= (1+2*lambda1);
        upp_diag(1:nx-2) = -lambda1;
        lower_diag(2:nx-1) = -lambda1;
        
        %defining the RHS matrix
        R = zeros(nx-1,1);
        R(1) = T_half(2,j) + lambda2*(T_half(2,j+1) - 2*(T_half(2,j)) + T_half(2,j-1))+ 75*lambda1;
        R(nx-1)= T_half(nx,j) + lambda2*(T_half(nx,j+1) - 2*(T_half(nx,j)) + T_half(nx,j-1))+ 50*lambda1;
        for i=2:nx-2
            R(i) = T_half(i+1,j) + lambda2*(T_half(i+1,j+1) - 2*(T_half(i+1,j)) + T_half(i+1,j-1));
        end
        %applying thomas algo
        X = thomas(upp_diag,lower_diag,main_diag,R,nx-1);
        %final values of temp2
        for k=2:nx
            Temp2(k,j,n+1)=X(k-1);
        end
    end
end

%graph plotting 
figure;
[M,N] = size(Temp2(:,:,4*t));
[x,y] = meshgrid(1:N,1:M); 

%surface plot
subplot(2,2,1);
surf(x,y,Temp2(:,:,4*t));
colorbar;
view(0,90);

%counter plot
subplot(2,2,3);
contour(x,y,Temp2(:,:,4*t));
colorbar;

%surface plot
subplot(2,2,2);
surf(x,y,Temp2(:,:,4*t));
colorbar;

%surface plot different view
subplot(2,2,4);
surf(x,y,Temp2(:,:,4*t));
view(0,15);
colorbar;
sgtitle("Using ADI Method")


%defining thomas function 
function X=thomas(A,B,D,R,n)
    A(1) = A(1)/D(1);
    R(1) = R(1)/D(1);
    for i=2:n-1
        A(i)= A(i)/(D(i)- (B(i)*A(i-1)));
        R(i)= (R(i)-(B(i)*R(i-1)))/(D(i)- (B(i)*A(i-1)));
    end
    R(n)=(R(n)-(B(n)*(R(n-1))))/(D(n)- (B(n)*A(n-1)));
    ans = zeros(n,1);
    ans(n)=R(n);
    for i=n-1:-1:1
        ans(i) = R(i) - A(i)*(ans(i+1));
    end
    X =ans;
end
