% Second implementation 
% Utilisation du schema de Newton-Rhafson

%% ---------------- Initialisation des paramètres ---------------- %%
% Chargement des constantes physiques
physical_constants;


% Définition des paramètres de simulation
simulation_parameters;

% Propriétés des matériaux et tension de bande plate
material_properties;

% Paramètres de la jonction PN
pn_parameters;


%% ---------------- Définition du maillage ---------------- %%
% Maillage global
[X, Nx] = maillage(Wtot, dx); % X est le vecteur des positions et Nx son nombre de points

% /!\ ATTENTION : Nx est mis à jour et peut être différent de la version initiale

% Initialisation des vecteurs de charge et de potentiel
n = zeros(size(X));  % Densité d'électrons
p = zeros(size(X));  % Densité de trous
Vx = zeros(size(X)); % Potentiel électrique

% Initialisation du profil de charge (hypothèse de jonction abrupte)
rho_init = charge_initialisation(X, Na, Nd, Wp, Wn);

Vg=0;
Vd=Vbi;
Vb= Vbi*ones(size(X));
% Résolution initiale de l'équation de Poisson
V_init = Poisson1(X, rho_init, Vg, Vd, Nx);


Vprec = V_init;
[n, p, rho] = charge_classique(V_init, Vb, X, Na, Nd);

Erreur =1; % Initialisation de l'erreur de convergence


% Boucle de convergence de la résolution de Poisson
iter=0;
while (Erreur > SetErr)
    % Calcul du Laplacien de Poisson
    Lp = Laplacien_Poisson(dx, Nx);
    [LpNR, RHS_NR] = Poisson_NR(Lp, n, p, rho, Vprec);


    % Application des conditions aux limites (tension appliquée)
    [LpNR, RHS_NR] = boundary_cond(LpNR, RHS_NR, dx, Vg, Vd); % ici on est a l'equilibre donc on ne prend que la Vbi pour les conditions aux limites

    % Résolution de l'équation de Poisson
    V_sol = LpNR \ RHS_NR'; 

    [n,p,rho] = charge_classique(V_sol,Vb,X,Na,Nd,ni);   % Mise à jour des profils de charge avec la methode classique

    % Calcul de l'erreur quadratique et mise à jour du potentiel
    
    clc;
    Erreur = rms(V_sol - Vprec); 

    Vprec = V_sol;

    % Affichage des résultats sous forme de graphiques
    % Plot Results
    figure;
    subplot(3,1,1);
    plot(X, V_sol, 'r', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Potentiel \psi(x) (V)');
    title('Electrostatic Potentiel');
    grid on;

    subplot(3,1,2);
    semilogy(X, n, 'b', X, p, 'g', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Carrier Densities (m^{-3})');
    legend('Electrons n(x)', 'Holes p(x)');
    title('Carrier Densities');
    grid on;

    subplot(3,1,3);
    plot(X, rho, 'k', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Charge Density \rho(x) (C/m^3)');
    title('Space Charge Density');
    grid on;

    pause(0.1)

    
    iter=iter+1;

end

toc