clear all
close all
clc

rho = 1;
Cp = 1;

Lx = 10;
Ly = 10;

Nx = 51;
Ny  = 51;
Nt = 500;
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);

K = ones(Nx,Ny);
K(30:35, 20:25) = 0.001;

Sx = round(7*Nx/Lx);
Sy = round(3*Ny/Ly);

c = 1;
C = 0.05;    %C*dx*Nt=T
dt = C*dx/c;

T = zeros(Nx,Ny);

%Boundary Conditions
T(:,1) = 0;  %bottom
T(: , end) = 0; %top
T(1,:) = 0  ;%left
T(end,:) = T(end-1,:);%neumann


%initial conditions
T(:,:) = 0;
t = 0;



for n=1:Nt
    
    Told = T;
    
    t = t+dt;
    
    
    for j = 2:Nx-1
        for i = 2:Ny-1
            T(i,j) = Told(i,j) + dt*(K(i,j)/(rho/Cp))*((Told(i+1,j) + Told(i,j+1) - 4*Told(i,j) + Told(i-1,j) + Told(i,j-1))/dx/dx);
        end
    end
    
    %introducing source term for 3 sec
    
    if (t<=3)
        T(Sx,Sy) = T(Sx,Sy) + dt*1000/(rho*Cp);
    end
    
    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);
    [X,Y] = meshgrid(x,y);
    mesh(x,y,T);
    axis([0 Lx 0 Ly 0 50]);
    title(sprintf('Time = %f seconds',t));
    pause(0.01);
    
end


