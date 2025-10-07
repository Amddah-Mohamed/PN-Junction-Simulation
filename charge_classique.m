function  [n,p,rho] = charge_classique(Vinc,Vb,X,Na,Nd,ni)

physical_constants;
simulation_parameters;
%%% A COMPLETER %%%
    % Indixation du cote P (left) et du cote N(right)
    left_indices = (X <= 0);
    right_indices = (X > 0);
    V=Vb-Vinc;
    % mise a jour des densites de porteurs 
    n(left_indices) = (ni^2/Na).* exp( Vinc(left_indices) ./ kT);
    p(left_indices) = Na .* exp(- Vinc(left_indices) ./ kT);
    
    n(right_indices) = Nd .* exp( -V(right_indices) ./ kT);
    p(right_indices) = (ni^2/Nd).* exp( V(right_indices) ./ kT);

    % mise a jour de la densite de charges 
    rho = q*(p - n + Nd .* right_indices - Na .* left_indices);
    
end