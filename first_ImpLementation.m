%Première implémentation de la résolution de l’équation de Poisson

clc
% Load Physical Constants
physical_constants;

% Define Simulation Parameters
simulation_parameters;

% Define Material Properties and Flat-band Voltage
material_properties;

% Define PN Junction Parameters
pn_parameters;

% Define the Grid
[X, Nx] = maillage(Wtot, dx);

% Initialize Charge Density
rho_init = charge_initialisation(X, Na, Nd, Wp, Wn);
rho = rho_init; % Ensure consistency

% Set Convergence Criteria
Limit = 1e-16;

% Initialize Carrier Densities
n = zeros(size(X));
p = zeros(size(X));

% Define Boundary Conditions
V1 = 0;
V2 = Vbi;
Vb= Vbi*ones(size(X));
Vprev=zeros(size(X));

iter=0;
while true % Infinite loop, exits when condition is met
    iter=iter+1;
    % Solve Poisson Equation for Vx
    Vx = Poisson1(X, rho, V1, V2, Nx);

    % Update Carrier Densities Using Boltzmann Approximation
    [n,p,rho_new] = charge_classique(Vx,Vb,X,Na,Nd,ni);

    % Plot Results
    figure;
    subplot(3,1,1);
    plot(X, Vx, 'r', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Potentiel \psi(x) (V)');
    title('Electrostatic Potentiel');
    grid on;

    subplot(3,1,2);
    plot(X, n, 'b', X, p, 'g', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Carrier Densities (m^{-3})');
    legend('Electrons n(x)', 'Holes p(x)');
    title('Carrier Densities');
    grid on;

    subplot(3,1,3);
    plot(X, rho_new, 'k', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Charge Density \rho(x) (C/m^3)');
    title('Space Charge Density');
    grid on;

    pause(0.1)
    
    % Convergence
    if rms(Vx - Vprev) < Limit
        break; % Stop the loop if converged or max iterations reached
    end
    Vprev=Vx;
    rho = rho_new; % Update charge density for next iteration
end

fprintf('Converge en %d iterations.\n', iter);
