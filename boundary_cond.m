function [LpNR, RHS_NR] = boundary_cond(LpNR, RHS_NR, dx, Vg, Vd)
    % Applique les conditions aux limites pour la jonction PN sous polarisation directe
   

    N = length(RHS_NR);  % Nombre de points du maillage

    % Appliquer la condition de Dirichlet à la borne gauche (p-side)
    
    RHS_NR(1) = RHS_NR(1)-Vg/dx^2;    

    % Appliquer la condition de Dirichlet à la borne droite (n-side)
   
    RHS_NR(N) = RHS_NR(N)-Vd/(dx^2);
    

end


