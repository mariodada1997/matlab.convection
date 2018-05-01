% Script: convect2d.m
close all;
clear all;
 
% Specify x range and number of points
x0 = -3;
x1 = 3;
Nx = 40;
 
% Specify y range and number of points
y0 = -3;
y1 = 3;
Ny = 40;
 
% Construct mesh
x = linspace(x0,x1,Nx+1);
y = linspace(y0,y1,Ny+1);
[xg,yg] = ndgrid(x,y);
 
% Construct mesh needed for plotting
xp = zeros(4,Nx*Ny);
yp = zeros(4,Nx*Ny);
n = 0;
for j = 1:Ny
for i = 1:Nx
n = n + 1;
xp(1,n) = x(i);
yp(1,n) = y(j);
xp(2,n) = x(i+1);
yp(2,n) = y(j);
xp(3,n) = x(i+1);
yp(3,n) = y(j+1);
xp(4,n) = x(i);
yp(4,n) = y(j+1);
end
end
 
% Calculate midpoint values in each control volume
xmid = 0.5*(x(1:Nx) + x(2:Nx+1));
ymid = 0.5*(y(1:Ny) + y(2:Ny+1));
[xmidg,ymidg] = ndgrid(xmid,ymid);
 
% Calculate cell size in control volumes (assumed equal) and area
dx = x(2) - x(1);
dy = y(2) - y(1);
A = dx*dy;
 
% Set velocity
u = 1;
v = 1;
 
% Set final time
tfinal = 1;
 
% Set timestep
CFL = 1.0;
dt = 0.075;%CFL/(abs(u)/dx + abs(v)/dy);
 
% Set initial condition to Q0 = exp(-x^2 - 20*y^2)
% Note: technically, we should average the initial
% distribution in each cell but I chose to just set
% the value of Q in each control volume to the midpoint
% value of Q0.
Q = exp(-xmidg.^2 - 20*ymidg.^2);
t = 0;
 
% Loop until t > tfinal
while (t < tfinal)
% The following implement the bcs by creating a larger array
% for Q and putting the appropriate values in the first and last
% columns or rows to set the correct bcs
Qbc(2:Nx+1,2:Ny+1) = Q; % Copy Q into Qbc
 
%Setting periodic bc in x-direction
Qbc( 1,2:Ny+1) = Q(Nx, :); % Periodic bc
Qbc(Nx+2,2:Ny+1) = Q( 1, :); % Periodic bc
 
%Set periodic bc in y-direction
Qbc( 2:Nx+1, 1) = Q(:, Ny);
Qbc(2:Nx+1, Ny+2) = Q( :, 1);
 
% Calculate the flux at each interface
% First the i interfaces
F_left = -0.5 *( Qbc(2:Nx+1, 2:Ny+1) + Qbc(1:Nx, 2:Ny+1)) + 0.5 *( Qbc(2:Nx+1, 2:Ny+1) - Qbc(1:Nx, 2:Ny+1));
F_right = 0.5 *( Qbc(2:Nx+1, 2:Ny+1) + Qbc(3:Nx+2, 2:Ny+1)) + 0.5 *( Qbc(2:Nx+1, 2:Ny+1) - Qbc(3:Nx+2, 2:Ny+1));
 
% Now the j interfaces
F_upper = 0.5 * ( Qbc(2:Nx+1, 2:Ny+1) + Qbc(2:Nx+1, 3:Ny+2)) + 0.5 * ( Qbc(2:Nx+1,2:Ny+1) - Qbc(2:Nx+1, 3:Ny+2));
F_lower = -0.5 * ( Qbc(2:Nx+1, 2:Ny+1) + Qbc(2:Nx+1, 1:Ny)) + 0.5 * ( Qbc(2:Nx+1,2:Ny+1) - Qbc(2:Nx+1, 1:Ny));
 
% Add contributions to right hand side from fluxes
F = F_right + F_left + F_upper + F_lower;
 
% Forward Euler step
Q = Q -(dt/dx) * F;
 
% Increment time
t = t + dt;
 
% Plot current solution
Qp = reshape(Q,1,Nx*Ny);
clf;
[Hp] = patch(xp,yp,Qp);
set(Hp,'EdgeAlpha',0);
axis('equal');
caxis([0,1]);
colorbar;
drawnow;
end
