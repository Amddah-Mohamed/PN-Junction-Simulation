function rho_init = charge_initialisation(X,Na,Nd,Wp,Wn)

    physical_constants
    
    rho_init=zeros(size(X));
    rho_init(X<=-Wp) = 0; % Région P neutre 
    rho_init(-Wp<X & X<=0) = -q*Na;  % Région P ZCE
    rho_init(0<X & X<=Wn) = q*Nd ;   %Region N ZCE
    rho_init(Wn<X) =0;   %Region N neutre 

end

