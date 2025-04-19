function [tout,yout] = RK4Solver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
% the right-hand side of the ODE ydot = f(t,y)
% t =[t0 t1 ... tfinal] is a vector of grid points
% with length N
% y0=[a; b; c] is a column vector that contain the
% initial values y(0) = y0.
% OUTPUT:tout is a column vector of grid points.
% yout is an 3 x N matrix containing the solution
% at different grid points.

% Initialization
tout = t(:);
m=length(y0);
N = length(tout);   
yout = zeros(m, N);
yout(:,1) = y0;

% Runge–Kutta (ERK) Methods
for n = 1:N-1

    tn = tout(n);         
    yn = yout(:,n);        
    h = tout(n+1) - tn; % find gap
    
    k1 = h * f(tn, yn);
    k2 = h * f(tn + h/2, yn + k1/2);
    k3 = h * f(tn + h/2, yn + k2/2);
    k4 = h * f(tn + h, yn + k3);
    yout(:,n+1) = yn + (k1 + 2*k2 + 2*k3 + k4)/6;
end
end