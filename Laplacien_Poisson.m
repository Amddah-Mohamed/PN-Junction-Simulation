function Lp = Laplacien_Poisson(dx,Nx)

material_properties;

% Calcul des diagonales supérieurs (Ap) et inférieur (Bp)

    % Question preliminaire : combien de point pour la matrice Laplacienne ?

%%% A COMPLETER %%%

  e = ones(Nx, 1);

  Lp = (1/dx^2)*spdiags([e -2*e e], -1:1, Nx, Nx);


end
