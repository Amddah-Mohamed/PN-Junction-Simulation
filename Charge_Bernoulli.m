% Mise a jour des densites a l'aide des fonctions de Bernoulli 
function [n, p, rho] = Charge_Bernoulli(V_sol, X,n ,p , T) %Le T et la duree de vie de porteur, pour ne pas le considerer il suffit de mettre T=0 

    % Importation des constantes et des paramètres du matériau
    material_properties;
    physical_constants;
    simulation_parameters;

    % Conversion du potentiel en unités normalisées par la tension thermique
    V = (1 / Vt) * V_sol;
    N = length(X);  % Nombre de points de grille

    % Initialisation des matrices et vecteurs du système
    An = zeros(N, N);  % Matrice du système pour les électrons
    Ap = zeros(N, N);  % Matrice du système pour les trous
    bn = zeros(N, 1);  % Second membre pour les électrons
    bp = zeros(N, 1);  % Second membre pour les trous
    
    % Construction du système par discrétisation 
    for i = 2:N-1
        % Construction de la matrice pour les électrons
        An(i, i-1) = B(V(i-1) - V(i));
        An(i, i)   =  -B(V(i) - V(i-1)) - B(V(i) - V(i+1));
        An(i, i+1) = B(-V(i) + V(i+1));

        % Construction de la matrice pour les trous
        Ap(i, i-1) = B(-V(i-1) + V(i));
        Ap(i, i)   = - B(-V(i) + V(i-1)) - B(-V(i) + V(i+1));
        Ap(i, i+1) = B(V(i) - V(i+1));
        
        if T>0 %Pour prendre la generation en consideration 
            bn(i)=(q*(n(i)*p(i)-ni^2))/(T*(n(i)+p(i)+2*ni));
            bp(i)=(-q*(n(i)*p(i)-ni^2))/(T*(n(i)+p(i)+2*ni));
        end 
    end

    % Normalisation des matrices par la mobilité et l'énergie thermique
    An = (q*muSi * kT / dx^2)* An;
    Ap = (-q*muSi * kT / dx^2)* Ap;

    % Application des conditions limites de Dirichlet (valeurs fixées aux bords)
    An(1, 1) = 1;
    An(N, N) = 1;
    Ap(1, 1) = 1;
    Ap(N, N) = 1;

    % Définition des valeurs aux frontières en fonction des concentrations de dopage
    bn(1) = (-Na + sqrt(Na^2 + 4 * ni^2)) / 2; % Frontière gauche (zone p)
    bn(N) = (Nd + sqrt(Na^2 + 4 * ni^2)) / 2;  % Frontière droite (zone n)
    bp(1) = (Na + sqrt(Na^2 + 4 * ni^2)) / 2;  % Frontière gauche (zone p)
    bp(N) = (-Nd + sqrt(Na^2 + 4 * ni^2)) / 2; % Frontière droite (zone n)

    % Résolution des systèmes pour obtenir n et p
    n = An \ bn;
    p = Ap \ bp;

    % Conversion des vecteurs en lignes pour correspondre au format souhaité
    n = n';
    p = p';

    % Détermination des indices correspondant aux régions dopées n et p
    left_indices = (X <= 0);   % Région p
    right_indices = (X > 0);   % Région n
    dopage = zeros(size(X));
    dopage(right_indices) = Nd;   % Dopage n à droite
    dopage(left_indices) = -Na;   % Dopage p à gauche

    % Mise à jour de la densité de charge rho
    rho = q * (p - n + dopage);
end
