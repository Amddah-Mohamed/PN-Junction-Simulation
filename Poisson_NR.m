function [LpNR,RHS_NR] = Poisson_NR(Lp,n,p,rho,Vprec)

    physical_constants;
    
    material_properties;
    
    simulation_parameters;
    
    % Chargement du Laplacien de Poisson
    
    LpNR=Lp-(q/(EpsSi*kT))*diag(n+p);

    
    
    % Surcharge de la diagonnale du Laplacien
    
    %%% A COMPLETER %%%
    
    % Modification du "Right Hand Side" (second membre) de l'équation de Poisson
    
       
    RHS_NR=(1/EpsSi)*(-rho-((1/kT)*Vprec').*(q*(n+p)));
        



end