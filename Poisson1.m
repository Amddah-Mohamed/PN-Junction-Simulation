function Vx = Poisson1(X,rho,V1,V2,Nx)
    

    material_properties ;

    dx=X(2)-X(1);
    % Création d'un vecteur d'éléments 1 pour les diagonales 
    %adjacentes
    e = ones(Nx, 1);
    
    % Construction de la matrice Laplacienne
    Laplacien = (1/dx^2)*spdiags([e -2*e e], -1:1, Nx, Nx);
    Y=-(rho/EpsSi);       
    % Conditions aux limites de Dirichlet
    Y(1) = Y(1)-V1/dx^2;
    
    % For last node (x=L)
    Y(Nx) =Y(Nx)- V2/dx^2;
    
    % Solve system
    Vx = Laplacien\Y';
end

