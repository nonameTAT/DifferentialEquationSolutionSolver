function [tout,yout] = EulerSolver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
% the right-hand side of the ODE ydot = f(t,y)
% t =[t0 t1 ... tfinal] is a vector of grid points
% with length N
% y0=[a; b; c] is a column vector that contain the
% initial values y(0)=y0.
% OUTPUT:tout is a column vector of grid points.
% yout is an 3 x N matrix containing the solution
% at different grid points.

% Initialization
tout = t(:);
N = length(tout);
m = length(y0);
yout = zeros(m, N);
yout(:,1) = y0;

% Euler's method
for n = 1:N-1
    tn = tout(n);
    yn = yout(:,n);
    h = tout(n+1) - tn;
    ydot = f(tn, yn);
    yout(:,n+1) = yn + h * ydot;
end
end