% 3eme implementation en Determinant les densites de charges avec l'aide
% des fonctions de Bernoulli 

clc; close all; clear all;
tic


% Chargement des constantes physiques
physical_constants;


% Définition des paramètres de simulation
simulation_parameters;

% Propriétés des matériaux et tension de bande plate
material_properties;

% Paramètres de la jonction PN
pn_parameters;


% Maillage global
[X, Nx] = maillage(Wtot, dx); % X est le vecteur des positions et Nx son nombre de points
% /!\ ATTENTION : Nx est mis à jour et peut être différent de la version initiale

% Initialisation des vecteurs de charge et de potentiel
n = zeros(size(X));  % Densité d'électrons
p = zeros(size(X));  % Densité de trous
Vx = zeros(size(X)); % Potentiel électrique

% Initialisation du profil de charge (hypothèse de jonction abrupte)
rho_init = charge_initialisation(X, Na, Nd, Wp, Wn);

% Conditions aux limites du potentiel
V1 = 0;       % Potentiel à gauche (polarisation appliquée)
V2 = Vbi;     % Potentiel à droite (tension de contact interne)
Vb = Vbi * ones(size(X));

% Résolution initiale de l'équation de Poisson
V_init = Poisson1(X, rho_init, V1, V2, Nx);



Erreur = 1; % Initialisation de l'erreur de convergence
Vprec = V_init;
[n, p, rho] = charge_classique(V_init, Vb, X, Na, Nd);

for j=1:length(Vapp_plus)


    Erreur =1;
    Vg=0;
    Vd=Vbi-Vapp_plus(j);

    % Boucle de convergence de la résolution de Poisson
    while (Erreur > SetErr)
        % Calcul du Laplacien de Poisson
        Lp = Laplacien_Poisson(dx, Nx);
        [LpNR, RHS_NR] = Poisson_NR(Lp, n, p, rho, Vprec);
        % Vg=0;
        % Vd=Vbi-Vapp(j);

        % Application des conditions aux limites (tension appliquée)
        [LpNR, RHS_NR] = boundary_cond(LpNR, RHS_NR, dx, Vg, Vd);

        % Résolution de l'équation de Poisson
        V_sol = LpNR \ RHS_NR';

        
        %[n,p,rho] = charge_classique(V_sol,Vd,X,Na,Nd,ni); % Mise à jour des profils de charge avec la methode classique 
        %SEulement pour montrer qu'elle ne donne pas des solutions pour J
                                                           

        [n, p, rho] = Charge_Bernoulli(V_sol, X,n,p,0);   % Mise à jour des profils de charge avec la fonction de Bernoulli

        % Calcul de l'erreur quadratique et mise à jour du potentiel
        clc;
        Erreur = rms(V_sol - Vprec); 

        Vprec = V_sol;

      

    end

    [XJ, J, Jn, Jp] = Courants(V_sol, X, n, p); %determination de J
    
    % Affichage des résultats sous forme de graphiques

    figure;
    subplot(4,1,1);
    plot(X, V_sol, 'r', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Potentiel \psi(x) (V)');
    title(sprintf('Electrostatic Potentiel for Va= %f', Vapp_plus(j)));
    grid on;

    subplot(4,1,2);
    semilogy(X, n, 'b', X, p, 'g', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Carrier Densities (m^{-3})');
    legend('Electrons n(x)', 'Holes p(x)');
    title('Carrier Densities');
    grid on;

    subplot(4,1,3);
    plot(X, rho, 'k', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Charge Density \rho(x) (C/m^3)');
    title('Space Charge Density');
    grid on;

    subplot(4,1,4);
    plot(XJ, J, 'g', 'LineWidth', 2);
    xlabel('Position x (m)');
    ylabel('Charge Density \rho(x) (C/m^3)');
    title('Space Charge Density');
    grid on;



    pause(0.1)

end
   

toc
