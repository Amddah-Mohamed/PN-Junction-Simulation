%% ---------------------------------------------------- %%
%% ---------------- Simulation de la jonction PN ------ %%
%% ---------------------------------------------------- %%
clc; close all; clear all;
tic

%% ---------------- Initialisation des param�tres ---------------- %%
% Chargement des constantes physiques
physical_constants;


% D�finition des param�tres de simulation
simulation_parameters;

% Propri�t�s des mat�riaux et tension de bande plate
material_properties;

% Param�tres de la jonction PN
pn_parameters;


%% ---------------- D�finition du maillage ---------------- %%
% Maillage global
[X, Nx] = maillage(Wtot, dx); % X est le vecteur des positions et Nx son nombre de points
% /!\ ATTENTION : Nx est mis � jour et peut �tre diff�rent de la version initiale

% Initialisation des vecteurs de charge et de potentiel
n = zeros(size(X));  % Densit� d'�lectrons
p = zeros(size(X));  % Densit� de trous
Vx = zeros(size(X)); % Potentiel �lectrique

% Initialisation du profil de charge (hypoth�se de jonction abrupte)
rho_init = charge_initialisation(X, Na, Nd, Wp, Wn);

% Conditions aux limites du potentiel
V1 = 0;       % Potentiel � gauche (polarisation appliqu�e)
V2 = Vbi;     % Potentiel � droite (tension de contact interne)
Vb = Vbi * ones(size(X));

% R�solution initiale de l'�quation de Poisson
V_init = Poisson1(X, rho_init, V1, V2, Nx);



I_plus=zeros(size(Vapp_plus));

I_moins=zeros(size(Vapp_moins));


%% ---------------- Boucle de convergence PN � l'�quilibre ---------------- %%

Erreur = 1; % Initialisation de l'erreur de convergence
Vprec = V_init;
[n, p, rho] = charge_classique(V_init, Vb, X, Na, Nd);


colors = lines(length(Tau_list)); % G�n�rer un ensemble de couleurs distinctes

for i=1:length(Tau_list)

    for j=1:length(Vapp_plus)


        Erreur =1;
        Vg=0;
        Vd=Vbi-Vapp_plus(j);

        % Boucle de convergence de la r�solution de Poisson
        while (Erreur > SetErr)
            % Calcul du Laplacien de Poisson
            Lp = Laplacien_Poisson(dx, Nx);
            [LpNR, RHS_NR] = Poisson_NR(Lp, n, p, rho, Vprec);
            % Vg=0;
            % Vd=Vbi-Vapp(j);

            % Application des conditions aux limites (tension appliqu�e)
            [LpNR, RHS_NR] = boundary_cond(LpNR, RHS_NR, dx, Vg, Vd);

            % R�solution de l'�quation de Poisson
            V_sol = LpNR \ RHS_NR'; %%% � COMPLETER %%%

            %[n,p,rho] = charge_classique(V_sol,Vb,X,Na,Nd,ni);   % Mise � jour des profils de charge avec la methode classique


            [n, p, rho] = Charge_Bernoulli(V_sol, X,n,p,Tau_list(i));   % Mise � jour des profils de charge avec la fonction de Bernoulli

            % Calcul de l'erreur quadratique et mise � jour du potentiel
            clc;
            Erreur = rms(V_sol - Vprec); %%% � COMPLETER %%%

            Vprec = V_sol;

        end

        [XJ, J, Jn, Jp] = Courants(V_sol, X, n, p);

        I_plus(j)=abs(mean(J));

    end


    Vprec = V_init;
    [n, p, rho] = charge_classique(V_init, Vb, X, Na, Nd);

    for j=1:length(Vapp_moins)


        Erreur =1;
        Vg=0;
        Vd=Vbi-Vapp_moins(length(Vapp_moins)-j+1);
        % Boucle de convergence de la r�solution de Poisson
        while (Erreur > SetErr)
            % Calcul du Laplacien de Poisson
            Lp = Laplacien_Poisson(dx, Nx);
            [LpNR, RHS_NR] = Poisson_NR(Lp, n, p, rho, Vprec);
            % Vg=0;
            % Vd=Vbi-Vapp(j);

            % Application des conditions aux limites (tension appliqu�e)
            [LpNR, RHS_NR] = boundary_cond(LpNR, RHS_NR, dx, Vg, Vd);

            % R�solution de l'�quation de Poisson
            V_sol = LpNR \ RHS_NR'; %%% � COMPLETER %%%

            %[n,p,rho] = charge_classique(V_sol,Vb,X,Na,Nd,ni);   % Mise � jour des profils de charge avec la methode classique


            [n, p, rho] = Charge_Bernoulli(V_sol, X,n,p,Tau_list(i));   % Mise � jour des profils de charge avec la fonction de Bernoulli

            % Calcul de l'erreur quadratique et mise � jour du potentiel
            clc;
            Erreur = rms(V_sol - Vprec); %%% � COMPLETER %%%

            Vprec = V_sol;

        end

        [XJ, J, Jn, Jp] = Courants(V_sol, X, n, p);

        I_moins(length(Vapp_moins)-j+1)=abs(mean(J));

    end

    I=[I_moins,I_plus];

    % Affichage des r�sultats sous forme de graphiques
    figure;

    semilogy(Vapp, I, 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', sprintf('\\tau = %.1e', Tau_list(i)));
    xlabel('Tension appliqu�e V_{app} (V)');
    ylabel('Courant I (A)');
    legend show; % Afficher la l�gende
    title('Caract�ristique I(V) pour diff�rentes valeurs de \Tau');
    grid on;

    
end

toc
