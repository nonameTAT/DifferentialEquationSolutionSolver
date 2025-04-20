function [tout,yout] = IRK4Solver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
% the right-hand side of the ODE ydot = f(t,y)
% t =[t0 t1 ... tfinal] is a vector of grid points
% with length N
% y0 is a column vector that contains the initial values y(0) = y0.
% OUTPUT: tout is a column vector of grid points.
% yout is an 3 x N matrix containing the solution at different grid points.

% Coefficients
sqrt3 = sqrt(3);
c1 = 0.5 - sqrt3/6;
c2 = 0.5 + sqrt3/6;
A11 = 1/4;
A12 = 1/4 - sqrt3/6;
A21 = 1/4 + sqrt3/6;
A22 = 1/4;
b1 = 1/2;
b2 = 1/2;

N = length(t);
yout = zeros(length(y0), N);
yout(:,1) = y0;
tout = t(:); % Ensure is a column vector

for n = 1:N-1
    current_t = t(n);
    current_y = yout(:,n);
    h = t(n+1) - current_t;
    
    % define the system F(K) = 0
    t1 = current_t + c1*h;
    t2 = current_t + c2*h;
    
    % initial guess for k1 and k2: evaluate f at t1 and t2 with current_y
    k1_initial = f(t1, current_y);
    k2_initial = f(t2, current_y);
    K_initial = [k1_initial; k2_initial];
    
    % Define the function F for fsolve
    F = @(K) [
        K(1:length(y0)) - f(t1, current_y + h*(A11*K(1:length(y0)) + A12*K(length(y0)+1:end)));...
        K(length(y0)+1:end) - f(t2, current_y + h*(A21*K(1:length(y0)) + A22*K(length(y0)+1:end)))
    ];
    
    % Solve function
    options = optimset('Display', 'off');
    K_sol = fsolve(F, K_initial, options);
    
    % Find k1 and k2
    k1 = K_sol(1:length(y0));
    k2 = K_sol(length(y0)+1:end);
    
    % Next y
    y_next = current_y + h*(b1*k1 + b2*k2);
    
    % Store the result
    yout(:,n+1) = y_next;
end
end