% Set values for the parameters ρ, σ, β.
rho = 28;
sigma = 10;
beta = 8/3;

% Define the right hand side of the Lorenz system (1) as an
% anonymous function. It should follow the form
f = @(t,x) lorenz(t,x,sigma,rho,beta);

k = input('Enter value of k: 1,2,3,4,5,6\n')
h = 10^(-k);
tfinal = 1;
t = [0:h:tfinal];

y0 = [-1;3;4];
    
% handle unavailable choice
while 1
    fprintf('Method 1: EulerSolver\n');
    fprintf('Method 2: RK4Solver\n');
    fprintf('Method 3: IRK4Solver\n');
    method = input('Number of method to solve system 1: \n')
    
    if method == 1 || method == 2 || method == 3
        break
    end
    fprintf('invalid number, enter again\n\n');
end

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
title('3D plot of System');
grid on;