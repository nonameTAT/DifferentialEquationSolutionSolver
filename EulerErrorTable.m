% Set values for the parameters ρ, σ, β.
rho = 28;
sigma = 10;
beta = 8/3;

% Define the right hand side of the Lorenz system as an
% anonymous function. It should follow the form
f = @(t,x) lorenz(t,x,sigma,rho,beta);

for k = 2:6
    h = 10^(-k);
    tfinal = 1;
    t = [0:h:tfinal];
    
    y0 = [-1;3;4];
            
    [~,Y] = EulerSolver(f,t,y0);
    
    options = odeset('RelTol', 3.1e-14, 'AbsTol', 1e-16);
    [~, Ym] = ode45(f, t, y0, options); % Ym is Nx3
    error = max(max(abs(Y - Ym')));
    fprintf('k=%d, max error is %d\n', k, error);
end
