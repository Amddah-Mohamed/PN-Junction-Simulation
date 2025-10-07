% Determination des courants et de la caracteristique

%% ---------------------------------------------------- %%
%% ---------------- Simulation de la jonction PN ------ %%
%% ---------------------------------------------------- %%
clc; close all; clear all;
tic

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

% Conditions aux limites du potentiel
V1 = 0;       % Potentiel à gauche (polarisation appliquée)
V2 = Vbi;     % Potentiel à droite (tension de contact interne)
Vb = Vbi * ones(size(X));

% Résolution initiale de l'équation de Poisson
V_init = Poisson1(X, rho_init, V1, V2, Nx);





I=zeros(size(Vapp_plus));



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
        V_sol = LpNR \ RHS_NR'; %%% À COMPLETER %%%

        %[n,p,rho] = charge_classique(V_sol,Vb,X,Na,Nd,ni);   % Mise à jour des profils de charge avec la methode classique


        [n, p, rho] = Charge_Bernoulli(V_sol, X,n,p,0);   % Mise à jour des profils de charge avec la fonction de Bernoulli

        % Calcul de l'erreur quadratique et mise à jour du potentiel
        clc;
        Erreur = rms(V_sol - Vprec); %%% À COMPLETER %%%

        Vprec = V_sol;

    end

    [XJ, J, Jn, Jp] = Courants(V_sol, X, n, p);

    I(j)=abs(mean(J));

end



% Affichage des résultats sous forme de graphiques
figure;

plot (Vapp_plus, I, 'g', 'LineWidth', 2);
xlabel('Tension appliquée V_{app} (V)');
ylabel('Courant I (A/m^{2})');
title('Caractéristique I(V) ');
grid on;

figure;

semilogy(Vapp_plus, I, 'g', 'LineWidth', 2);
xlabel('Tension appliquée V_{app} (V)');
ylabel('Courant I (A)');
title('Caractéristique I(V) en echelle log ');
grid on;


toc
