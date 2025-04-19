% Set values for the parameters ρ, σ, β.
rho = 14;
sigma = 10;
beta = 8/3;

% Define the right hand side of the Lorenz system (1) as an
% anonymous function. It should follow the form
f = @(t,x) lorenz(t,x,sigma,rho,beta);

k = 3;
h = 10^(-k);
tfinal = 100;
t = [0:h:tfinal];

y0 = [-1;3;4];
    
method = 2;

if method == 1
    [~,Y] = EulerSolver(f,t,y0);
elseif method == 2
    [~,Y] = RK4Solver(f,t,y0);
elseif method == 3
    [~,Y] = IRK4Solver(f,t,y0);
end

options = odeset('RelTol', 3.1e-14, 'AbsTol', 1e-16);
[~, Ym] = ode45(f, t, y0, options); % Ym is Nx3
error = max(max(abs(Y - Ym')));
disp(error);

figure;
plot3(Y(1, :), Y(2, :), Y(3, :));
xlabel('x');
ylabel('y');
zlabel('z');
title('3D Trajectory of the Lorenz System');
grid on;